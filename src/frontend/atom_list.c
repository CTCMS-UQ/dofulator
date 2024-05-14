#include <stdlib.h>

#include "atom_list.h"

AtomList atom_list_create(unsigned n) {
  AtomList atoms = {.n = n};
  atoms.x = malloc(sizeof((*atoms.x)) * n);
  if (atoms.x) {
    atoms.mass = malloc(sizeof((*atoms.mass)) * n);
    if (atoms.mass) {
      return atoms;
    }
  }
  return (AtomList){.n = 0};
}

void atom_list_destroy(AtomList* atoms) {
  if (atoms->x) free(atoms->x);
  if (atoms->mass) free(atoms->mass);
  *atoms = (AtomList){.n = 0};
}
