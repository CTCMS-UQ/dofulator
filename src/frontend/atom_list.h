#ifndef DOF_ATOM_LIST_H
#define DOF_ATOM_LIST_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * A list of point masses
*/
typedef struct AtomList {
  unsigned n;      // Number of atoms
  double (*x)[3];  // Positions of each atom
  double *mass;    // Masses of each atom
} AtomList;

/*
 * Allocate an atom list with space for `n` atoms.
*/
AtomList atom_list_create(unsigned n);

/*
 * Clean up an atom list
*/
void atom_list_destroy(AtomList* atoms);

#ifdef __cplusplus
}
#endif

#endif
