from libc.stddef import size_t

cdef extern from "dofulator.h":
    ctypedef void* Dofulator
    ctypedef size_t AtomTag

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

    Dofulator dofulator_create(AtomTag n_atoms)
    void dofulator_destroy(Dofulator* ctx)

    void dofulator_add_rigid_bond(Dofulator ctx, Bond b)
    void dofulator_build_rigid_fragment(Dofulator ctx, Bond b)
    void dofulator_finalise_fragments(Dofulator ctx)

    void dofulator_precalculate_rigid(Dofulator ctx, const double* mass, const double x[][3])
    void dofulator_calculate(Dofulator ctx, const double* mass, const double x[][3])

    void dofulator_get_dof_atom_directional(const void* ctx, AtomTag atom_idx, double dof[3])
    double dofulator_get_dof_atom(const void* ctx, AtomTag atom_idx)

    FragmentListIter dofulator_get_rigid_fragments(const void* ctx);
    FragmentListIter dofulator_get_semirigid_fragments(const void* ctx);

    const Fragment* fragmentlist_iter_next(FragmentListIter* iter);

    AtomListView fragment_get_atoms(const Fragment* frag);
