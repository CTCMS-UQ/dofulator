#ifndef LIBDOFULATOR_H
#define LIBDOFULATOR_H

#include "atom_list.h"
#include "bond_list.h"

#ifdef __cplusplus
extern "C" {
#endif

// A rigid or semi-rigid fragment of `Atom`s connected by constraints
typedef void* Fragment;

/*
 * Create a semi-rigid fragment from an atom list and bond list.
 * Treats all bonds as fixed length with flexible angles, and assumes all atoms
 * are connected in a single graph.
 *
 * Bond list will be sanitized (mutated)
 *
 * WARNING: Currently cannot handle rings - only tree topology!
*/
Fragment fragment_create_semirigid(AtomList atoms, BondList bonds);

/*
 * Create a rigid fragment from an atom list.
 * All atoms treated as part of a single rigid body.
*/
Fragment fragment_create_rigid(AtomList atoms);

/*
 * Clean up a fragment.
*/
void fragment_destroy(Fragment fragment);

/*
 * Calculate the total DoF of the fragment
*/
double fragment_dof(const Fragment fragment);

/*
 * Calculate the DoF of atom with index `atom`
*/
double fragment_dof_atom(Fragment fragment, unsigned atom);

/*
 * Calculate the DoF of atom with index `atom` in direction `dir`.
 * `dir` must be a unit vector.
*/
double fragment_dof_atom_dir(Fragment fragment, unsigned atom, const double dir[3]);

#ifdef __cplusplus
}
#endif

#endif
