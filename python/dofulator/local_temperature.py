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
        verbose: bool = True,
        store_dof_results: bool = False,
    ):
        """
        MDAnalysis extension for calculating temperatures of groups of atoms.

        `selections` contains a list of AtomGroups, each of which will have
        its temperature calculated on each frame.

        Uses MDADofulator to calculate DoF in each group.
        AtomGroup used to construct `dofulator` should cover all atoms which
        may be in any of the groups in `selections`.
        If the direct results from `dofulator` are also required, set
        `store_dof_results = True`. Otherwise, they will be discarded after
        each frame.
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
            self._dof_buf = np.zeros((max((np.max(sel.ix) for sel in self.selections if len(sel) > 0)),))
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
                sel.masses * np.sum(sel.velocities * sel.velocities, axis=1)
            ) / dof

    def _conclude(self):
        self.results /= self.boltz

