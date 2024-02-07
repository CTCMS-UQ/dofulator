#ifndef DOF_BOND_LIST_H
#define DOF_BOND_LIST_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * A bond between atoms with indices i and j
*/
typedef struct Bond {
  unsigned i;
  unsigned j;
} Bond;

/*
 * List of bonds
*/
typedef struct BondList {
  unsigned n;   // Number of bonds
  Bond* bonds;  // List of atom index pairs
} BondList;

/*
 * Allocate memory for a bond list with `n` bonds
*/
BondList bond_list_create(unsigned n);

/*
 * Clean up a bond list
*/
void bond_list_destroy(BondList* bonds);

#ifdef __cplusplus
}
#endif

#endif
