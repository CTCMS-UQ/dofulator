import numpy as np
from typing import Iterable

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase

from .mda_dofulator import MDADofulator

class LocalTemperature(AnalysisBase):
    def __init__(
        self,
        selections: Iterable[mda.AtomGroup],
        dof: MDADofulator|np.ndarray|None = None,
        mode: str = 'atomic',
        dof_modes: str = 'all',
        store_dof_results: bool = False,
        verbose: bool = True,
    ):
        """
        MDAnalysis extension for calculating temperatures of groups of atoms.
        Assumes kinetic energy of each atom is commensurate with constraints.

        After running, the `results` member will be a numpy array with `n_frames` rows
        and `len(selections)` columns, where each column holds the temperature
        on each frame of the selection with index equal to the column index.

        Note: If MDAnalysis incorrectly guesses velocity units, the `boltz` member may
        need to be modified.

        Parameters
        ----------
        `selections`: `Iterable[AtomGroup]`
            A list of AtomGroups, each of which will have its temperature
            calculated on each frame. For spatial selections, make sure to set
            `updating=True` in the `AtomGroup`.

        `dof`: `None|NDArray|MDADofulator`
            `None`: All atoms assigned 3 DoF.
            `NDArray`: Element `ix` should have the total DoF of atom with index `ix`.
            `MDADofulator`: should be initialised to calculate DoF forc all atoms which
                could be in any of the groups in `selections`.

        `mode`: `'atomic'|'directional'`
            Passed through to the `MDADofulator` instance.
            'atomic': Only use atomic DoF and calculate just the total kinetic temperature
                 of each selection.
            'directional': Use directional DoF to calculate x,y,z kinetic temperatures.

        `dof_modes'`: `'all'|'translational'|'rovibrational'`
            Passed through to `MDADofulator` instance - determines which modes to include
            for temperature calculation.
            'all': Use all modes
            'translational': Use translational modes only (temperature calculated from CoM velocity)
            'rovibrational': Use rotational + vibrational modes only (temperature calculated from v - v_com)

        `store_dof_results`: `bool` (default False)
            Set `True` if the direct results from `dofulator` are also required,
            otherwise they will be discarded after each frame.
        """
        if type(selections) is list:
            self.selections = selections
        else:
            self.selections = list(selections)
        super(LocalTemperature, self).__init__(
            self.selections[0].universe.trajectory,
            verbose=verbose
        )
        self.dof = dof
        self.mode = mode
        self.dof_modes = dof_modes
        self.store_dof_results = store_dof_results
        self.boltz = mda.units.constants['Boltzmann_constant']

    def _prepare(self):
        if isinstance(self.dof, MDADofulator):
            self.dof.n_frames = self.n_frames if self.store_dof_results else 1
            self.dof.mode = self.mode
            self.dof._prepare()
        if self.mode == 'atomic':
            self.results = np.zeros((self.n_frames, len(self.selections)), dtype=np.float64)
            self._single_frame = self._single_frame_atomic
        elif self.mode == 'directional':
            self.results = np.zeros((self.n_frames, len(self.selections), 3), dtype=np.float64)
            self._single_frame = self._single_frame_directional
        else:
            raise Exception(f"Invalid mode '{self.mode}'. Must be one of 'atomic' or 'directional'")

    def _single_frame_atomic(self):
        from .c_dofulator import modes_from_str
        if isinstance(self.dof, MDADofulator):
            self._dof_buf = np.zeros((max((sel.n_atoms for sel in self.selections)),))
            if self.store_dof_results:
                self.dof._frame_index = self._frame_index
                self.dof._single_frame()
            else:
                # Just call _calculate if not storing results to avoid
                # duplicate work querying atomic DoF
                self.dof._calculate()

        if self.dof_modes != 'all':
            # Remove translational or rotational/vibrational component from velocities
            # Translational DoF is 3x mass fraction, so use that to sum up CoM velocity
            # TODO: this could be made faster by moving to C/Cython
            if not isinstance(self.dof, MDADofulator):
                raise Exception("dof_modes must be 'all' if DoF is not calculated by an `MDADofulator` instance")
            dof_trans = self.dof._ctx.get_all_dof(mode=modes_from_str('translational')) / 3
            velocities = self.dof._atomgroup.velocities.copy()
            if self.dof_modes == 'translational':
                for frag in self.dof._ctx.get_fragments():
                    v = np.zeros((3,))
                    for a in frag.atoms():
                        v += dof_trans[a] * velocities[a,:]
                    for a in frag.atoms():
                        velocities[a,:] = v
            else:
                for frag in self.dof._ctx.get_fragments():
                    v = np.zeros((3,))
                    for a in frag.atoms():
                        v += dof_trans[a] * velocities[a,:]
                    for a in frag.atoms():
                        velocities[a,:] -= v

        for i, sel in enumerate(self.selections):
            if self.dof is None:
                # No DoF provided. Assume 3 per atom.
                dof = 3 * len(sel)
            elif isinstance(self.dof, MDADofulator):
                # DoF from dofulator
                dof = np.sum(self.dof._ctx.get_all_dof(sel.ix, self._dof_buf, modes_from_str(self.dof_modes)))
            else:
                # DoF from array of constants
                dof = np.sum(self.dof[sel.ix])

            # Calculate m*v*v / d. Apply Boltzmann factor later.
            if self.dof_modes == 'all':
                self.results[self._frame_index, i] = np.sum(
                    sel.masses * np.sum(sel.velocities**2, axis=1)
                ) / dof
            else:
                self.results[self._frame_index, i] = np.sum(
                    sel.masses * np.sum(velocities[sel.ix]**2, axis=1)
                ) / dof


    def _single_frame_directional(self):
        from .c_dofulator import modes_from_str
        if isinstance(self.dof, MDADofulator):
            self._dof_buf = np.zeros((max((sel.n_atoms for sel in self.selections)),3))
            if self.store_dof_results:
                self.dof._frame_index = self._frame_index
                self.dof._single_frame()
            else:
                # Just call _calculate if not storing results to avoid
                # duplicate work querying atomic DoF
                self.dof._calculate()

        if self.dof_modes != 'all':
            # Remove translational or rotational/vibrational component from velocities
            # Translational DoF is 3x mass fraction, so use that to sum up CoM velocity
            # TODO: this could be made faster by moving to C/Cython
            if not isinstance(self.dof, MDADofulator):
                raise Exception("dof_modes must be 'all' if DoF is not calculated by an `MDADofulator` instance")
            dof_trans = self.dof._ctx.get_all_dof(mode=modes_from_str('translational')) / 3
            velocities = self.dof._atomgroup.velocities.copy()
            if self.dof_modes == 'translational':
                for frag in self.dof._ctx.get_fragments():
                    v = np.zeros((3,))
                    for a in frag.atoms():
                        v += dof_trans[a] * velocities[a,:]
                    for a in frag.atoms():
                        velocities[a,:] = v
            else:
                for frag in self.dof._ctx.get_fragments():
                    v = np.zeros((3,))
                    for a in frag.atoms():
                        v += dof_trans[a] * velocities[a,:]
                    for a in frag.atoms():
                        velocities[a,:] -= v

        for i, sel in enumerate(self.selections):
            if self.dof is None:
                # No DoF provided. Assume 3 per atom.
                dof = [len(sel)] * 3
            elif isinstance(self.dof, MDADofulator):
                # DoF from dofulator
                dof = np.sum(self.dof._ctx.get_all_dof_directional(sel.ix, self._dof_buf, modes_from_str(self.dof_modes)), axis=0)
            else:
                # DoF from array of constants not supported since values will quickly vary
                raise Exception("Directional temperature with constant DoF not supported. "
                                "Molecular rotation and deformation will quickly change the x,y,z "
                                "contributions, so they should be recalculated each step.")

            # Calculate m*v*v / d. Apply Boltzmann factor later.
            mass = np.repeat(np.reshape(sel.masses, (len(sel.masses), 1)), 3, axis=1)
            if self.dof_modes == 'all':
                self.results[self._frame_index, i] = np.sum(mass * sel.velocities**2, axis=0) / dof
            else:
                self.results[self._frame_index, i] = np.sum(mass * velocities[sel.ix]**2, axis=0) / dof

    def _conclude(self):
        self.results /= self.boltz

