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
        self.store_dof_results = store_dof_results
        self.boltz = mda.units.constants['Boltzmann_constant']

    def _prepare(self):
        if isinstance(self.dof, MDADofulator):
            self.dof.n_frames = self.n_frames if self.store_dof_results else 1
            self.dof._prepare()
        self.results = np.zeros((self.n_frames, len(self.selections)), dtype=np.float64)

    def _single_frame(self):
        if isinstance(self.dof, MDADofulator):
            self._dof_buf = np.zeros((max((sel.n_atoms for sel in self.selections)),))
            if self.store_dof_results:
                self.dof._frame_index = self._frame_index
                self.dof._single_frame()
            else:
                # Just call _calculate if not storing results to avoid
                # duplicate work querying atomic DoF
                self.dof._calculate()

        for i, sel in enumerate(self.selections):
            if self.dof is None:
                # No DoF provided. Assume 3 per atom.
                dof = 3 * len(sel)
            elif isinstance(self.dof, MDADofulator):
                # DoF from dofulator
                dof = np.sum(self.dof._ctx.get_all_dof(sel.ix, self._dof_buf))
            else:
                # DoF from array of constants
                dof = np.sum(self.dof[sel.ix])
            # Calculate m*v*v / d. Apply Boltzmann factor later.
            self.results[self._frame_index, i] = np.sum(
                sel.masses * np.sum(sel.velocities**2, axis=1)
            ) / dof

    def _conclude(self):
        self.results /= self.boltz

