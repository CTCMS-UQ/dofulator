import numpy as np
from typing import Iterable

import MDAnalysis as mda
from MDAnalysis.lib import mdamath
from MDAnalysis.core.groups import AtomGroup, Residue, ResidueGroup
from MDAnalysis.core.topologyobjects import Angle, Bond, Dihedral
from MDAnalysis.core.universe import Atom
from MDAnalysis.analysis.base import AnalysisBase

from .c_dofulator import CDofulator as Dofulator



class MDADofulator(AnalysisBase):
    def __init__(
        self,
        atomgroup: AtomGroup,
        rigid_bodies: str|ResidueGroup|Iterable[Iterable[Atom]]|None = None,
        rigid_bonds: str|Iterable[Bond]|None = None,
        rigid_angles: str|Iterable[Angle]|None = None,
        rigid_dihedrals: str|Iterable[Dihedral]|None = None,
        mode: str = 'atomic',
        verbose: bool = True,
        null_space_thresh: float|None = None,
        use_pbc: bool = True,
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

        `rigid_bodies`, `rigid_bonds`, `rigid_angles`, `rigid_dihedrals`: Topology groups,
        lists of topology objects, or `'all'`
        Sets features to treat as rigid. See `set_rigid_bodies()`, `set_rigid_bonds()`,
        `set_rigid_angles()` and `set_rigid_dihedrals()` for details.

        `mode`: `'atomic'|'directional'` sets whether to calculate only the
        total DoF per atom (`atomic`) or the DoF in each Cartesian direction
        ('directional')

        `null_space_thresh`: Threshold for calculating null space of loop closure matrix.
        Sanitised to be between 0. and 1. using `min(abs(thresh), 1.)`.
        `None|0` = use default from core C library (DBL_EPSILON * largest dimension of loop closure matrix).

        `use_pbc`: `bool` set `False` to ignore periodic boundaries.
        """
        super(MDADofulator, self).__init__(atomgroup.universe.trajectory, verbose=verbose)
        self._atomgroup = atomgroup
        self.mode = mode
        self.null_space_thresh = null_space_thresh
        self.use_pbc = use_pbc
        self.set_rigid_bodies(rigid_bodies)
        self.set_rigid_bonds(rigid_bonds)
        self.set_rigid_angles(rigid_angles)
        self.set_rigid_dihedrals(rigid_dihedrals)

    @property
    def null_space_thresh(self):
        if self._null_space_thresh is not None:
            return self._null_space_thresh
        elif hasattr(self, '_ctx') and self._ctx:
            return self._ctx.null_space_thresh
        else:
            # No context constructed yet, so use a dummy one to get the default
            return Dofulator(0).null_space_thresh

    @null_space_thresh.setter
    def null_space_thresh(self, thresh: float|None):
        self._null_space_thresh = min(abs(thresh), 1.) if thresh is not None else None


    def set_rigid_bodies(self, rigid_bodies: str|ResidueGroup|Iterable[Residue|Iterable[Atom]]|None = None):
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
        if isinstance(rigid_bodies, str):
            if rigid_bodies == 'all':
                self._rigid_bodies = [r.atoms for r in self._atomgroup.residues]
            else:
                raise ValueError(f'Invalid value "{rigid_bodies}" for argument `rigid_bodies`')
        elif isinstance(rigid_bodies, list) and len(rigid_bodies) > 0 \
        and (isinstance(rigid_bodies[0], list) or isinstance(rigid_bodies[0], AtomGroup)):
            # List of lists or list of AtomGroups - avoid deep copy
            # Assuming uniform list.
            self._rigid_bodies = rigid_bodies
        else:
            # Generator expression or something else unknown
            #   - check for Residues, or otherwise assume an atom list
            self._rigid_bodies = [
                    [a for a in (r.atoms if isinstance(r, Residue) else r)]
                    for r in rigid_bodies
            ] if rigid_bodies else None

    def set_rigid_bonds(self, rigid_bonds: str|Iterable[Bond]|None = None):
        """
        Set the list of rigid bonds.
        Atoms connected to each other by a set of rigid bonds form a semi-rigid fragment.

        `rigid_bonds`:
            * `'all'` treats every bond in the atomgroup as rigid.
            * Anything else is iterated over, with each item expected to be of type `Bond`.
        """
        if isinstance(rigid_bonds, str):
            if rigid_bonds == 'all':
                self._rigid_bonds = self._atomgroup.bonds
            else:
                raise ValueError(f'Invalid value "{rigid_bonds}" for argument `rigid_bonds`')
        elif isinstance(rigid_bonds, list):
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
        if isinstance(rigid_angles, str):
            if rigid_angles == 'all':
                self._rigid_angles = self._atomgroup.angles
            else:
                raise ValueError(f'Invalid value "{rigid_angles}" for argument `rigid_angles`')
        elif isinstance(rigid_angles, list):
            # Avoid deep copy of list
            self._rigid_angles = rigid_angles
        else:
            self._rigid_angles = [a for a in rigid_angles] if rigid_angles else None

    def set_rigid_dihedrals(self, rigid_dihedrals: str|Iterable[Dihedral]|None = None):
        """
        Set the list of rigid dihedrals.
        For simplicity, it is assumed (and automatically enforced) that the
        bonds and angles within the dihedral are also rigid. This behaviour may be
        relaxed in future versions.
        """
        if isinstance(rigid_dihedrals, str):
            if rigid_dihedrals == 'all':
                self._rigid_dihedrals = self._atomgroup.dihedrals
            else:
                raise ValueError(f'Invalid value "{rigid_dihedrals}" for argument `rigid_dihedrals`')
        elif isinstance(rigid_dihedrals, list):
            # Avoid deep copy of list
            self._rigid_dihedrals = rigid_dihedrals
        else:
            self._rigid_dihedrals = [d for d in rigid_dihedrals] if rigid_dihedrals else None

    def _setup_ctx(self):
        """
        Set up the Dofulator context for the current set of rigid bodies, bonds and angles.
        Any previously used context is destroyed.
        """
        max_ix = np.max(self._atomgroup.universe.atoms.ix)
        if hasattr(self, '_ctx') and self._ctx:
            del self._ctx
        self._ctx = Dofulator(max_ix+1)
        if self._null_space_thresh is not None:
            self._ctx.null_space_thresh = self._null_space_thresh

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
                # Construct bonds first for consistency
                self._ctx.add_rigid_bond(a.atoms[0].ix, a.atoms[1].ix)
                self._ctx.add_rigid_bond(a.atoms[1].ix, a.atoms[2].ix)
                self._ctx.add_rigid_bond(a.atoms[0].ix, a.atoms[2].ix)

        if self._rigid_dihedrals:
            for d in self._rigid_dihedrals:
                # Construct bonds first, then angles, then dihedral for consistency
                self._ctx.add_rigid_bond(d.atoms[0].ix, d.atoms[1].ix)
                self._ctx.add_rigid_bond(d.atoms[1].ix, d.atoms[2].ix)
                self._ctx.add_rigid_bond(d.atoms[2].ix, d.atoms[3].ix)
                self._ctx.add_rigid_bond(d.atoms[0].ix, d.atoms[2].ix)
                self._ctx.add_rigid_bond(d.atoms[1].ix, d.atoms[3].ix)
                self._ctx.add_rigid_bond(d.atoms[0].ix, d.atoms[3].ix)

        self._ctx.finalise_fragments()

    def _set_ctx_pbc(self):
        """
        Update PBCs in dofulator context.
        """
        if not self.use_pbc:
            self._ctx.set_pbc_none()
            return

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
        else:
            raise Exception(f"Invalid mode '{self.mode}'. Must be one of 'atomic' or 'directional'")
        all_atoms = self._atomgroup.universe.atoms
        self._masses = all_atoms.masses.astype(np.float64)  # Assume masses don't change during run
        self._set_ctx_pbc()
        self._ctx.precalculate_rigid(self._masses, all_atoms.positions.astype(np.float64))

    def _calculate(self):
        # PBCs may be changing, so update them
        self._set_ctx_pbc()
        self._ctx.calculate(self._masses, self._atomgroup.universe.atoms.positions.astype(np.float64))

    def _single_frame(self):
        self._calculate()
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

