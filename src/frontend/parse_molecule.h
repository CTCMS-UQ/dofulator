#ifndef DOF_PARSE_MOLECULES_H
#define DOF_PARSE_MOLECULES_H

#include "atom_list.h"
#include "bond_list.h"

typedef struct Molecule {
  char (*atom_labels)[8];
  AtomList atoms;
  BondList bonds;
} Molecule;

Molecule parse_molecule_file(char* fname);

void molecule_destroy(Molecule* mol);

#endif
