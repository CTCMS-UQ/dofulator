import numpy as np
from typing import Iterable

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase

from .mda_dofulator import MDADofulator

class LocalTemperature(AnalysisBase):
    def __init__(
        self,
        selections: Iterable[mda.AtomGroup],
        dofulator: MDADofulator|None = None,
        store_dof_results: bool = False,
        verbose: bool = True,
    ):
        """
        MDAnalysis extension for calculating temperatures of groups of atoms.
        Assumes kinetic energy of each atom is commensurate with constraints.

        `selections`: a list of AtomGroups, each of which will have
        its temperature calculated on each frame. For spatial selections,
        make sure to set `updating=True` in the `AtomGroup`.

        `dofulator`: Can be a list of DoF values for all atoms in the
        `Universe`, or an `MDADofulator`, in which case it should cover
        all atoms which could be in any of the groups in `selections`.

        `store_dof_results`: (default False) Set `True` if the direct results
        from `dofulator` are also required otherwise, they will be discarded
        after each frame.

        Note: If MDAnalysis incorrectly guesses velocity units, the `boltz` member may
        need to be modified.
        """
        if type(selections) is list:
            self.selections = selections
        else:
            self.selections = list(selections)
        super(LocalTemperature, self).__init__(
            self.selections[0].universe.trajectory,
            verbose=verbose
        )
        self.dofulator = dofulator
        self.store_dof_results = store_dof_results
        self.boltz = mda.units.constants['Boltzmann_constant']

    def _prepare(self):
        if self.dofulator is not None:
            self.dofulator.n_frames = self.n_frames if self.store_dof_results else 1
            self.dofulator._prepare()
        self.results = np.zeros((self.n_frames, len(self.selections)), dtype=np.float64)

    def _single_frame(self):
        if self.dofulator is not None:
            self._dof_buf = np.zeros((max((sel.n_atoms for sel in self.selections)),))
            if self.store_dof_results:
                self.dofulator._frame_index = self._frame_index
                self.dofulator._single_frame()
            else:
                # Just call _calculate if not storing results to avoid
                # duplicate work querying atomic DoF
                self.dofulator._calculate()

        for i, sel in enumerate(self.selections):
            if self.dofulator is None:
                dof = 3 * len(sel)
            else:
                dof = np.sum(self.dofulator._ctx.get_all_dof(sel.ix, self._dof_buf))
            # Calculate m*v*v / d. Apply Boltzmann factor later.
            self.results[self._frame_index, i] = np.sum(
                sel.masses * np.sum(sel.velocities**2, axis=1)
            ) / dof

    def _conclude(self):
        self.results /= self.boltz

