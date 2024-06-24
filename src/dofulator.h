#ifndef LIBDOFULATOR_H
#define LIBDOFULATOR_H

#include <stdlib.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// May be overridden by build system for compatibility
#ifndef AtomTag
#define AtomTag int64_t
#endif

// TODO: This should be left at 1 for now since batch operations
//       aren't implemented properly yet.
#ifndef DOFULATOR_DEFAULT_BATCH_SIZE
#define DOFULATOR_DEFAULT_BATCH_SIZE 1
#endif

/*******************************************************************************
 * Opaque handle to dofulator context
*/
typedef struct Dofulator* Dofulator;

/*******************************************************************************
 * One fragment. Public API only deals with pointers to this
*/
typedef struct Fragment Fragment;

/*******************************************************************************
 * Bond between atoms with indices `ai` and `aj`
*/
typedef struct Bond {
  AtomTag ai;
  AtomTag aj;
} Bond;

/*******************************************************************************
 * Iterator over a list of fragments.
 * May no longer be valid if underlying fragment list is mutated.
*/
typedef struct FragmentListIter {
  size_t n_fragments;         // Number of fragments to iterate over
  size_t idx;                 // Current index
  const Fragment* fragments;  // Read-only view of fragments. Invalidated if FragmentList changes.
} FragmentListIter;

/*******************************************************************************
 * Read-only view of a list of atoms.
 * May no longer be valid if underlying atom list is mutated.
*/
typedef struct AtomListView {
  size_t n_atoms;         // Number of atoms.
  const AtomTag* atoms;   // Atom indices.
                          // Read-only, but may be changed/invalidated if
                          // fragments are rebuilt. Should be copied out
                          // if persistence is needed.
} AtomListView;


/*******************************************************************************
 * Information about simulation box/periodic boundaries.
*/
typedef struct PBC {
  enum {
    PBC_NONE = 0,   // Non-periodic box
    PBC_TRI,        // Triclinic box. a along x axis, b in xy plane.
    PBC_ORTHO       // Orthogonal box.
  } typ;
  union {
    struct {
      double lx, _pad_yzx[3], ly, _pad_zxy[3], lz;
    };
    struct {
      double ax, _pad_a[2];
      double bx, by, _pad_b;
      double cx, cy, cz;
    };
    struct {
      double a[3], b[3], c[3];
    };
  };
} PBC;


/*******************************************************************************
 * Context management
*/

// Create a dofulator context with space for `n_atoms` atoms
Dofulator dofulator_create(AtomTag n_atoms);

// Clean up a dofulator context.
void dofulator_destroy(Dofulator* ctx);

// Set current boundary conditions
void dofulator_set_pbc(Dofulator ctx, PBC pbc);



/*******************************************************************************
 * Preparation
*/

// Add a rigid bond. This will create/modify/merge semi-rigid fragments as needed.
// Must call `dofulator_finalise_fragments(ctx)` once all fragments are built.
void dofulator_add_rigid_bond(Dofulator ctx, Bond b);

// Mark the two atoms in the bond as being in the same rigid body.
// Will add/modify/merge rigid fragments as needed.
// Must call `dofulator_finalise_fragments(ctx)` once all fragments are built.
void dofulator_build_rigid_fragment(Dofulator ctx, Bond b);

// Finalise fragment construction and allocate working memory.
// No fragments should be added or modified after this point.
void dofulator_finalise_fragments(Dofulator ctx);



/*******************************************************************************
 * Calculation
*/

// Pre-calculate DoF and reference frame for rigid fragments.
// Assumes mass will be constant through future calculations, and that the geometry in `x`
// is representative of the rigid bodies (i.e. constraints are respected).
// May be re-called with updated masses if necessary.
void dofulator_precalculate_rigid(Dofulator ctx, const double* mass, const double x[][3]);

// Calculate directional DoF of all atoms in semi-rigid fragments, and relative
// orientations of rigid bodies (from which directional DoF can be obtained later if needed)
// Assumes `dofulator_precalculate_rigid()` has been called.
void dofulator_calculate(Dofulator ctx, const double* mass, const double x[][3]);

// Get the total DoF of atom with index `atom_idx`.
// Assumes `dofulator_calculate()` has been called.
double dofulator_get_dof_atom(const struct Dofulator* ctx, AtomTag atom_idx);

// Get the directional DoF of atom with index `atom_idx` (returned in `dof` as `{x,y,z}`).
// Assumes `dofulator_calculate()` has been called.
void dofulator_get_dof_atom_directional(const struct Dofulator* ctx, AtomTag atom_idx, double dof[3]);

// Get the total DoF of the `n_atoms` atoms with indices in `atoms`. Returned in `dof`.
// Assumes `dofulator_calculate()` has been called.
double* dofulator_get_dof_atoms(
  const struct Dofulator* ctx, const size_t n_atoms, const AtomTag* atoms, double* dof);

// Get the directional DoF of the `n_atoms` atoms with indices in `atoms`. Returned in `dof`.
// Assumes `dofulator_calculate()` has been called.
double* dofulator_get_dof_atoms_directional(
  const struct Dofulator* ctx, const size_t n_atoms, const AtomTag* atoms, double dof[][3]);



/*******************************************************************************
 * Inspection
*/

// Get an iterator over the list of rigid bodies
FragmentListIter dofulator_get_rigid_fragments(const struct Dofulator* ctx);

// Get an iterator over the list of semi-rigid fragments
FragmentListIter dofulator_get_semirigid_fragments(const struct Dofulator* ctx);

// Get the next fragment from the iterator
const Fragment* fragmentlist_iter_next(FragmentListIter* iter);

// Get a read-only view of atoms in a given fragment
AtomListView fragment_get_atoms(const Fragment* frag);



#ifdef __cplusplus
}
#endif

#endif
