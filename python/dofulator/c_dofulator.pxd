from libc.stddef import size_t
from libc.stdint import int64_t

cdef extern from "dofulator.h":
    ctypedef void* Dofulator
    ctypedef signed long AtomTag

    cdef struct Bond:
        AtomTag ai
        AtomTag aj

    cdef struct Fragment:
        pass
    
    cdef struct FragmentListIter:
        size_t n_fragments
        size_t idx
        const Fragment* fragments

    cdef struct AtomListView:
        size_t n_atoms
        const AtomTag* atoms

    cdef enum PBC_typ:
        PBC_NONE
        PBC_TRI
        PBC_ORTHO
    cdef struct PBC:
        PBC_typ typ
        # Below are packed into an anonymous union,
        # but Cython doesn't know about memory packing
        # so this works to access them
        double lx, ly, lz
        double ax, bx, by, cx, cy, cz
        double a[3]
        double b[3]
        double c[3]
        # To prevent warnings:
        double _pad_a[2]
        double _pad_b
        double _pad_yzx[3]
        double _pad_zxy[3]

    Dofulator dofulator_create(AtomTag n_atoms)
    void dofulator_destroy(Dofulator* ctx)
    void dofulator_set_pbc(Dofulator ctx, PBC pbc)

    void dofulator_add_rigid_bond(Dofulator ctx, Bond b)
    void dofulator_build_rigid_fragment(Dofulator ctx, Bond b)
    void dofulator_finalise_fragments(Dofulator ctx)

    void dofulator_precalculate_rigid(Dofulator ctx, const double* mass, const double x[][3])
    void dofulator_calculate(Dofulator ctx, const double* mass, const double x[][3])

    void dofulator_get_dof_atom_directional(const void* ctx, AtomTag atom_idx, double dof[3])
    double dofulator_get_dof_atom(const void* ctx, AtomTag atom_idx)

    double* dofulator_get_dof_atoms(void* ctx, const size_t n_atoms, const AtomTag* atoms, double* dof);
    double* dofulator_get_dof_atoms_directional(void* ctx, const size_t n_atoms, const AtomTag* atoms, double dof[][3]);

    FragmentListIter dofulator_get_rigid_fragments(const void* ctx);
    FragmentListIter dofulator_get_semirigid_fragments(const void* ctx);

    const Fragment* fragmentlist_iter_next(FragmentListIter* iter);

    AtomListView fragment_get_atoms(const Fragment* frag);
