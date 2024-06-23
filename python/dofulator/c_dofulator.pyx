cimport cython
from cython cimport view

import numpy
cimport numpy
numpy.import_array()

cimport dofulator.c_dofulator as dof
from dofulator.c_dofulator cimport Dofulator, AtomTag, Bond, FragmentListIter, AtomListView

cdef class CDofulator:
    cdef Dofulator ctx
    cdef AtomTag n_atoms
    def __cinit__(self, AtomTag n_atoms):
        self.ctx = dof.dofulator_create(n_atoms)
        self.n_atoms = n_atoms

    def __dealloc__(self):
        dof.dofulator_destroy(&self.ctx)

    def add_rigid_bond(self, AtomTag i, AtomTag j):
        """
        Add a rigid bond constraint between atoms i and j

        NOTE: make sure `finalise_fragments()` is called once
        all fragments have been built.
        """
        if i >= self.n_atoms or j >= self.n_atoms:
            raise IndexError("Atom index out of range")
        dof.dofulator_add_rigid_bond(self.ctx, Bond(i, j))

    def build_rigid_fragment(self, AtomTag i, AtomTag j):
        """
        Mark atoms i and j as being in the same rigid body

        NOTE: make sure `finalise_fragments()` is called once
        all fragments have been built.
        """
        if i >= self.n_atoms or j >= self.n_atoms:
            raise IndexError("Atom index out of range")
        dof.dofulator_build_rigid_fragment(self.ctx, Bond(i, j))

    def finalise_fragments(self):
        """
        Finalise all rigid and semi-rigid fragments ready for
        DoF calculation.
        """
        # TODO: Handle C memory errors
        dof.dofulator_finalise_fragments(self.ctx)

    def set_pbc_none(self):
        """
        Turn off periodic boundary conditions
        """
        cdef PBC pbc
        pbc.typ=PBC_NONE
        dofulator_set_pbc(self.ctx, pbc)

    def set_pbc_ortho(self, numpy.ndarray[numpy.float64_t,ndim=1] box):
        """
        Set periodic boundaries with an orthorhombix box.
        `box = [lx, ly, lz]`
        """
        cdef PBC pbc
        pbc.typ = PBC_ORTHO
        pbc.lx = box[0]
        pbc.ly = box[1]
        pbc.lz = box[2]
        dofulator_set_pbc(self.ctx, pbc)

    def set_pbc_triclinic(self, numpy.ndarray[numpy.float64_t,ndim=2] box):
        """
        Set periodic boundaries with a triclinic box.
        `box = [[ax, 0, 0], [bx, by, 0], [cx, cy, cz]]`
        """
        cdef PBC pbc
        pbc.typ = PBC_TRI
        pbc.ax = box[0,0]
        pbc.bx = box[1,0]
        pbc.by = box[1,1]
        pbc.cx = box[2,0]
        pbc.cy = box[2,1]
        pbc.cz = box[2,2]
        dofulator_set_pbc(self.ctx, pbc)

    def precalculate_rigid(
        self,
        numpy.ndarray[numpy.float64_t,ndim=1] mass,
        numpy.ndarray[numpy.float64_t,ndim=2,mode='c'] x
    ):
        """
        Calculate DoF of all rigid bodies. These are invariant, and can therefore
        be rotated with each body as needed.
        """
        if len(mass) < self.n_atoms:
            raise IndexError("mass array too short")
        if x.shape[0] < self.n_atoms:
            raise IndexError("x array too short")
        if x.shape[1] != 3:
            raise IndexError("x array must be n_atoms x 3")
        dof.dofulator_precalculate_rigid(self.ctx,
            cython.cast(cython.pointer(double), mass.data),
            cython.cast(cython.pointer(double[3]), x.data))

    def calculate(
        self,
        numpy.ndarray[numpy.float64_t,ndim=1] mass,
        numpy.ndarray[numpy.float64_t,ndim=2,mode='c'] x
    ):
        """
        Calculate DoF of all atoms in all semi-rigid fragments
        in the given conformation and relative rotations of
        all rigid fragments.
        Assumes `self.precalculate_rigid(...)` has been called.
        """
        if len(mass) < self.n_atoms:
            raise IndexError("mass array too short")
        if x.shape[0] < self.n_atoms:
            raise IndexError("x array too short")
        if x.shape[1] != 3:
            raise IndexError("x array must be n_atoms x 3")
        dof.dofulator_calculate(self.ctx,
            cython.cast(cython.pointer(double), mass.data),
            cython.cast(cython.pointer(double[3]), x.data))

    def get_dof(self, AtomTag i):
        """
        Get the total DoF of atom with index `i`.
        """
        if i >= self.n_atoms:
            raise IndexError
        return dof.dofulator_get_dof_atom(self.ctx, i)

    def get_dof_directional(self, AtomTag i):
        """
        Get the DoF of atom with index `i` in each Cartesian direction
        """
        if i >= self.n_atoms:
            raise IndexError
        cdef double[3] d
        dof.dofulator_get_dof_atom_directional(self.ctx, i, d)
        return (d[0], d[1], d[2])

    def get_rigid_fragments(self):
        """
        Get iterator over rigid fragments
        """
        frags = FragmentIter()
        frags._iter = dofulator_get_rigid_fragments(self.ctx)
        return frags

    def get_semirigid_fragments(self):
        """
        Get iterator over semirigid fragments
        """
        frags = FragmentIter()
        frags._iter = dofulator_get_semirigid_fragments(self.ctx)
        return frags

    def get_all_dof(self, atoms=None, buffer=None):
        """
        Get DoF of all atoms in `atoms`.
        Optionally store result directly in `buffer`, which should be of type `np.float64`
        If `atoms` is `None`, gets DoF of all atoms (index 0 to self.n_atoms-1)
        """
        cdef AtomTag N
        if atoms is None:
            N = self.n_atoms
        else:
            N = len(atoms)
        if buffer is None:
            buf = numpy.empty((N,), dtype=numpy.float64, order='C')
            self._get_all_dof(atoms, buf)
            return buf
        if buffer.shape[0] < N:
            raise IndexError("`buffer` must have space for at least as many items in `atoms`")
        self._get_all_dof(atoms, buffer)
        return buffer

    def _get_all_dof(self, numpy.ndarray[numpy.int64_t,ndim=1] atoms, numpy.ndarray[numpy.float64_t,ndim=1] buffer):
        cdef AtomTag i
        cdef AtomTag N
        cdef double* buf = cython.cast(cython.pointer(double), buffer.data)
        if atoms is None:
            dof.dofulator_get_dof_atoms(self.ctx, buffer.shape[0], cython.NULL, buf)
        else:
            dof.dofulator_get_dof_atoms(
                self.ctx,
                atoms.shape[0],
                cython.cast(cython.pointer(AtomTag), atoms.data),
                buf
            )

    def get_all_dof_directional(self, atoms=None, buffer=None):
        """
        Get directional DoF of all atoms in `atoms`.
        Optionally store result directly in `buffer`, which should be an n x 3 array of type `np.float64`
        If `atoms` is `None`, gets DoF of all atoms (index 0 to self.n_atoms-1)
        """
        cdef AtomTag N
        if atoms is None:
            N = self.n_atoms
        else:
            N = len(atoms)
        if buffer is None:
            buf = numpy.empty((N, 3), dtype=numpy.float64, order='C')
            self._get_all_dof_directional(atoms, buf)
            return buf
        if buffer.shape[1] != 3 or buffer.shape[0] < N:
            raise IndexError("`buffer` must have at least as many rows as `atoms`, and exactly 3 columns")
        self._get_all_dof_directional(atoms, buffer)
        return buffer

    def _get_all_dof_directional(self, numpy.ndarray[numpy.int64_t,ndim=1] atoms, numpy.ndarray[numpy.float64_t,ndim=2,mode='c'] buffer):
        cdef AtomTag i
        cdef AtomTag N
        cdef double[3]* buf = cython.cast(cython.pointer(double[3]), buffer.data)
        if atoms is None:
            dof.dofulator_get_dof_atoms_directional(self.ctx, buffer.shape[0], cython.NULL, buf)
        else:
            dof.dofulator_get_dof_atoms_directional(
                self.ctx,
                atoms.shape[0],
                cython.cast(cython.pointer(AtomTag), atoms.data),
                buf
            )


cdef class FragmentHandle:
    """
    Opaque handle to a fragment
    """
    cdef const Fragment* h
    def atoms(self):
        """
        Get a read-only view of atom indices in this fragment.
        View is invalidated if the fragment list is mutated (e.g. by
        `CDofulator.add_rigid_bond()`, `CDofulator.build_rigid_fragment()`,
        or `CDofulator.finalise_fragments()`)
        """
        cdef AtomListView atoms = fragment_get_atoms(self.h)
        cdef const AtomTag [::1] atoms_v = <AtomTag[:atoms.n_atoms:1]>atoms.atoms
        return atoms_v


cdef class FragmentIter:
    """
    Iterator over a fragment list.
    Becomes invalid if the fragment list is mutated (e.g. by
    `CDofulator.add_rigid_bond()`, `CDofulator.build_rigid_fragment()`,
    or `CDofulator.finalise_fragments()`)
    """
    cdef FragmentListIter _iter
    def __len__(self):
        return self._iter.n_fragments - self._iter.idx
    def __iter__(self):
        return self
    def __next__(self):
        f = FragmentHandle()
        f.h = fragmentlist_iter_next(&self._iter)
        if f.h:
            return f
        else:
            raise StopIteration

