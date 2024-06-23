import numpy as np
from typing import Iterable

import MDAnalysis as mda
from MDAnalysis.lib import mdamath
from MDAnalysis.core.groups import ResidueGroup
from MDAnalysis.core.topologyobjects import Angle, Bond
from MDAnalysis.core.universe import Atom
from MDAnalysis.analysis.base import AnalysisBase

from .c_dofulator import CDofulator as Dofulator



class MDADofulator(AnalysisBase):
    def __init__(
        self,
        atomgroup: mda.AtomGroup,
        rigid_bodies: str|ResidueGroup|Iterable[Iterable[Atom]]|None = None,
        rigid_bonds: str|Iterable[Bond]|None = None,
        rigid_angles: str|Iterable[Angle]|None = None,
        mode: str = 'atomic',
        verbose: bool = True,
    ):
        """
        MDAnalysis extension for calculating per-atom degrees of freedom given
        a set of rigid constraints.
        Groups of atoms within `rigid_bodies` are treated as individual rigid bodies.
        Atoms connected by bonds in `rigid_bonds` are treated as a semi-rigid fragment.
        Atoms forming an angle in `rigid_angles` also contribute to semi-rigid fragments.
        Currently, rigid angles are treated by setting both bonds as rigid, and
        adding a virtual 3rd bond to form a rigid triangle. In future, this may
        be relaxed to allow stretchable bonds with a rigid angle.

        `rigid_bodies`, `rigid_bonds`, `rigid_angles`: see `set_rigid_bodies()`, `set_rigid_bonds()`
        and `set_rigid_angles()` for details

        `mode`: `'atomic'|'directional'` sets whether to calculate only the
        total DoF per atom (`atomic`) or the DoF in each Cartesian direction
        ('directional')
        """
        super(MDADofulator, self).__init__(atomgroup.universe.trajectory, verbose=verbose)
        self._atomgroup = atomgroup
        self.mode = mode
        self.pbc: bool = atomgroup.universe.dimensions is not None
        self.set_rigid_bodies(rigid_bodies)
        self.set_rigid_bonds(rigid_bonds)
        self.set_rigid_angles(rigid_angles)

    def set_rigid_bodies(self, rigid_bodies: str|ResidueGroup|Iterable[Iterable[Atom]]|None = None):
        """
        Set the list of rigid bodies.
        Rigid bodies must not be connected to semi-rigid fragments by any rigid constraints.

        `rigid_bodies`:
            * `'all'` treats each residue in the atomgroup as a rigid body.
            * `ResidueGroup` treats each residue in the `ResidueGroup` as a rigid body.
            * Anything else is iterated over, with each item iterated as a list of atoms
              within a unique rigid body. Any atom that appears in multiple rigid bodies will
              result in those bodies being joined.
        """
        if type(rigid_bodies) is str:
            if rigid_bodies == 'all':
                self._rigid_bodies = [[a for a in r.atoms] for r in self._atomgroup.residues]
            else:
                raise ValueError(f'Invalid value "{rigid_bodies}" for argument `rigid_bodies`')
        elif type(rigid_bodies) is ResidueGroup:
            self._rigid_bodies = [[a for a in r.atoms] for r in rigid_bodies]
        elif type(rigid_bodies) is list and (type(rigid_bodies[0] is list) or type(rigid_bodies[0]) is mda.AtomGroup):
            # Avoid deep copy of list
            self._rigid_bodies = rigid_bodies
        else:
            self._rigid_bodies = [[a for a in r] for r in rigid_bodies] if rigid_bodies else None

    def set_rigid_bonds(self, rigid_bonds: str|Iterable[Bond]|None = None):
        """
        Set the list of rigid bonds.
        Atoms connected to each other by a set of rigid bonds form a semi-rigid fragment.

        `rigid_bonds`:
            * `'all'` treats every bond in the atomgroup as rigid.
            * Anything else is iterated over, with each item expected to be of type `Bond`.
        """
        if type(rigid_bonds) is str:
            if rigid_bonds == 'all':
                self._rigid_bonds = self._atomgroup.bonds
            else:
                raise ValueError(f'Invalid value "{rigid_bonds}" for argument `rigid_bonds`')
        elif type(rigid_bonds) is list:
            # Avoid deep copy of list
            self._rigid_bonds = rigid_bonds
        else:
            self._rigid_bonds = [b for b in rigid_bonds] if rigid_bonds else None

    def set_rigid_angles(self, rigid_angles: str|Iterable[Angle]|None = None):
        """
        Set the list of rigid angles.
        For simplicity, it is assumed (and automatically enforced) that the
        bonds forming each angle are also rigid. This behaviour may be
        relaxed in future versions.
        """
        if type(rigid_angles) is str:
            if rigid_angles == 'all':
                self._rigid_angles = self._atomgroup.angles
            else:
                raise ValueError(f'Invalid value "{rigid_angles}" for argument `rigid_angles`')
        elif type(rigid_angles) is list:
            # Avoid deep copy of list
            self._rigid_angles = rigid_angles
        else:
            self._rigid_angles = [a for a in rigid_angles] if rigid_angles else None

    def _setup_ctx(self):
        """
        Set up the Dofulator context for the current set of rigid bodies, bonds and angles.
        Any previously used context is destroyed.
        """
        max_ix = max((ix for ix in self._atomgroup.universe.atoms.ix))
        if hasattr(self, '_ctx') and self._ctx:
            del self._ctx
        self._ctx = Dofulator(max_ix+1)

        if self._rigid_bodies:
            for b in self._rigid_bodies:
                a0 = None
                for a in b:
                    if a0 is not None:
                        self._ctx.build_rigid_fragment(a0, a.ix)
                    else:
                        a0 = a.ix

        if self._rigid_bonds:
            for b in self._rigid_bonds:
                self._ctx.add_rigid_bond(b.atoms[0].ix, b.atoms[1].ix)

        if self._rigid_angles:
            for a in self._rigid_angles:
                # If any of these bonds have already been added, nothing
                # will change for that bond, so safe to just add all.
                self._ctx.add_rigid_bond(a.atoms[0].ix, a.atoms[1].ix)
                self._ctx.add_rigid_bond(a.atoms[0].ix, a.atoms[2].ix)
                self._ctx.add_rigid_bond(a.atoms[1].ix, a.atoms[2].ix)

        self._ctx.finalise_fragments()

    def _set_ctx_pbc(self):
        """
        Update PBCs in dofulator context.
        """
        dimensions = self._atomgroup.universe.dimensions
        if dimensions is None:
            mda.warnings.warn('No dimension data detected. Assuming non-periodic '
                          'system. Dofulator will produce incorrect results '
                          'if this system should be periodic.')
            self._ctx.set_pbc_none()
            return
        if np.any(dimensions[3:6] != 90.):
            self._ctx.set_pbc_ortho(dimensions.astype(np.float64))
        else:
            self._ctx.set_pbc_triclinic(mdamath.triclinic_vectors(dimensions).astype(np.float64))


    def _prepare(self):
        self._setup_ctx()
        if self.mode == 'atomic':
            self.results = np.zeros((self.n_frames, self._atomgroup.n_atoms), dtype=np.float64)
        elif self.mode == 'directional':
            self.results = np.zeros((self.n_frames, self._atomgroup.n_atoms, 3), dtype=np.float64)
        all_atoms = self._atomgroup.universe.atoms
        self._masses = all_atoms.masses.astype(np.float64)
        self._set_ctx_pbc()
        self._ctx.precalculate_rigid(self._masses, all_atoms.positions.astype(np.float64))

    def _single_frame(self):
        self._set_ctx_pbc()
        self._ctx.calculate(self._masses, self._atomgroup.universe.atoms.positions.astype(np.float64))
        if self.results.ndim == 2:
            # mode == 'atomic'
            self._ctx.get_all_dof(self._atomgroup.ix, self.results[self._frame_index, :])
        else:
            # mode == 'directional'
            self._ctx.get_all_dof_directional(self._atomgroup.ix,self.results[self._frame_index, :, :])

    def _conclude(self):
        pass

    def run_single_frame(self):
        self.n_frames = 1
        self._frame_index = 0
        self._prepare()
        self._single_frame()
        self._conclude()

