#include "dofulator.h"
#include <stdio.h>
#include <math.h>

#define PI 3.14159265378979
#define DEG (PI/180.0)

int main(int argc, char* argv[]) {
  AtomList atoms = {
    .n = 3,
    .mass = (double*)(double[]){15.999, 1.008, 1.008},
    .pos = (double*)(double[]){
      0.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      cos(109.47*DEG), sin(109.47*DEG), 0.0,
    },
  };
  Fragment frag = fragment_create_rigid(atoms);
  if (!frag) {
    fprintf(stderr, "Failed to create rigid body!");
    return 1;
  }
  double dof = 0.0;
  for (unsigned i = 0; i < atoms.n; ++i) {
    double d = fragment_dof_atom(frag, i);
    dof += d;
    printf("Atom %d:\t%g\n", i, d);
  }
  printf("Total:\t%g\n", dof);

  printf("\nz direction\n");
  dof = 0.0;
  for (unsigned i = 0; i < atoms.n; ++i) {
    double d = fragment_dof_atom_dir(frag, i, (double[3]){1.0, 0.0, 0.0});
    dof += d;
    printf("Atom %d:\t%g\n", i, d);
  }
  printf("Total:\t%g\n", dof);
  return 0;
}
