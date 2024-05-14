#include "bond_list.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

BondList bond_list_create(unsigned n) {
  Bond* bonds = malloc(sizeof(*bonds) * n);
  if (!bonds) {
    n = 0;
  }
  return (BondList){
    .n = n,
    .bonds = bonds,
  };
}

void bond_list_destroy(BondList* bonds) {
  if (bonds->bonds) free(bonds->bonds);
  *bonds = (BondList){.n = 0};
}

