#include <stdlib.h>

#include "atom_list.h"

AtomList atom_list_create(unsigned n) {
  AtomList atoms = {.n = n};
  atoms.pos = malloc(sizeof(double)*4*n);
  if (atoms.pos) {
    atoms.mass = &atoms.pos[3*n];
    return atoms;
  } else {
    return (AtomList){.n = 0};
  }
}

void atom_list_destroy(AtomList* atoms) {
  if (atoms->pos) free(atoms->pos);
  *atoms = (AtomList){.n = 0};
}
