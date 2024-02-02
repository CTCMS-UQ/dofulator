#ifndef LIBDOFULATOR_H
#define LIBDOFULATOR_H

// A list of point masses, formatted x0,y0,z0,x1,y1,z1,...
typedef struct AtomList {
  unsigned n;
  double *pos;
  double *mass;
} AtomList;

// A bond between atoms with indices i and j
typedef struct Bond {
  unsigned i;
  unsigned j;
} Bond;

typedef struct BondList {
  unsigned n;
  Bond* bonds;
} BondList;

// A rigid or semi-rigid fragment of `Atom`s connected by constraints
typedef void* Fragment;

// Create a semi-rigid fragment from an atom list and a connectivity list.
// Assumes all atoms are connected by one or more constraints.
// Assumes all listed bonds have a rigid length.
Fragment fragment_create_semirigid(AtomList atoms, BondList bonds);

// Create a rigid fragment from an atom list
Fragment fragment_create_rigid(AtomList atoms);

// Free resources associated with the fragment
void fragment_destroy(Fragment fragment);

// Get total DoF of the fragment
double fragment_dof(const Fragment fragment);

// Get total DoF of an atom (by index)
double fragment_dof_atom(Fragment fragment, unsigned atom);

// Get DoF of an atom (by index) in the provided direction.
// `dir` should be a unit vector.
double fragment_dof_atom_dir(Fragment fragment, unsigned atom, double dir[3]);

#endif
