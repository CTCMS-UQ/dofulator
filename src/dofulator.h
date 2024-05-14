#ifndef LIBDOFULATOR_H
#define LIBDOFULATOR_H

// #include "atom_list.h"
// #include "bond_list.h"
// #include "angle_list.h"

#include <stdlib.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// A rigid or semi-rigid fragment of `Atom`s connected by constraints
// typedef void* Fragment;

#ifndef AtomTag
#define AtomTag size_t
#endif

typedef struct Bond {
  AtomTag ai;
  AtomTag aj;
} Bond;

// TODO: This should be left at 1 for now since batch operations
//       aren't implemented properly yet.
#ifndef DOFULATOR_DEFAULT_BATCH_SIZE 
#define DOFULATOR_DEFAULT_BATCH_SIZE 1
#endif

typedef struct Dofulator* Dofulator;

Dofulator dofulator_create(AtomTag n_atoms);
void dofulator_destroy(Dofulator* ctx);
void dofulator_add_rigid_bond(Dofulator ctx, Bond b);
void dofulator_build_rigid_fragment(Dofulator ctx, Bond b);
void dofulator_finalise_fragments(Dofulator ctx);
void dofulator_precalculate_rigid(Dofulator ctx, const double* mass, const double (*x)[3]);
void dofulator_calculate(Dofulator ctx, const double* mass, const double (*x)[3]);
void dofulator_get_dof_atom_directional(const struct Dofulator* ctx, AtomTag atom_idx, double dof[3]);
double dofulator_get_dof_atom(const struct Dofulator* ctx, AtomTag atom_idx);

/*
 * Create a semi-rigid fragment from an atom list and bond list.
 * Treats all bonds as fixed length with flexible angles, and assumes all atoms
 * are connected in a single graph.
 *
 * Bond list will be sanitized (mutated)
 *
 * WARNING: Currently cannot handle rings - only tree topology!
*/
// Fragment fragment_create_semirigid(AtomList atoms, BondList bonds, AngleList angles);

/*
 * Create a rigid fragment from an atom list.
 * All atoms treated as part of a single rigid body.
*/
// Fragment fragment_create_rigid(AtomList atoms);

/*
 * Clean up a fragment.
*/
// void fragment_destroy(Fragment fragment);

/*
 * Calculate the total DoF of the fragment
*/
// double fragment_dof(const Fragment fragment);

/*
 * Calculate the DoF of atom with index `atom`
*/
// double fragment_dof_atom(Fragment fragment, unsigned atom);

/*
 * Calculate the DoF of atom with index `atom` in direction `dir`.
 * `dir` must be a unit vector.
*/
// double fragment_dof_atom_dir(Fragment fragment, unsigned atom, const double dir[3]);

#ifdef __cplusplus
}
#endif

#endif
