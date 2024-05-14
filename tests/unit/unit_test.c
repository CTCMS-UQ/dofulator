#include <stdio.h>
#include <stdlib.h>

#include "unit_test.h"

static unsigned digits(unsigned n) {
  unsigned d = 1;
  while (n >= 10) {
    n /= 10;
    ++d;
  }
  return d;
}

bool run_test(const Test* test) {
  bool result = true;
  Dofulator ctx = dofulator_create(test->atoms.n);
  if (!ctx) {
    printf("FAILED\n\tError creating dofulator context!\n");
    return false;
  }
  for (Bond* b = test->bonds.bonds; b < test->bonds.bonds + test->bonds.n; ++b) {
    switch (test->test_type) {
      case RIGID:
        dofulator_build_rigid_fragment(ctx, *b);
        break;
      case FLEX:
        dofulator_add_rigid_bond(ctx, *b);
    }
  }
  dofulator_finalise_fragments(ctx);
  dofulator_precalculate_rigid(ctx, test->atoms.mass, test->atoms.x);
  dofulator_calculate(ctx, test->atoms.mass, test->atoms.x);

  // double dof_total = 0.0;
  for (size_t i = 0; i < test->atoms.n; ++i) {
    double dof_atom = 0.0;
    double dof[3];
    dofulator_get_dof_atom_directional(ctx, i, dof);
    for (unsigned d = 0; d < 3; ++d) {
      if (!feql(dof[d], test->dof[i][d])) {
        if (result) printf("FAILED\n");
        printf("\t! Atom %ld expected %.16g DoF in %c, got %.16g\n",
               i, test->dof[i][d], (char[]){'x', 'y', 'z'}[d], dof[d]);
        result = false;
      }
      dof_atom += test->dof[i][d];
      // dof_total += test->dof[i][d];
    }
    double dof_atom_actual = dofulator_get_dof_atom(ctx, i);
    if (!feql(dof_atom, dof_atom_actual)) {
      if (result) printf("FAILED\n");
      printf("\t! Atom %ld expected %.16g DoF total, got %.16g\n",
             i, dof_atom, dof_atom_actual);
      result = false;
    }
  }
  // double dof = fragment_dof(frag);
  // if (!feql(dof_total, dof)) {
  //   if (result) printf("FAILED\n");
  //   printf("\t! Expected %.16g DoF total, got %.16g\n",
  //          dof_total, dof);
  //   result = false;
  // }

  dofulator_destroy(&ctx);
  return result;
}

int run_testset(const TestSet testset) {
  printf("\nRunning test %s...\n", testset.name);
  unsigned npassed = 0;
  const unsigned ndig = digits(testset.n);
  for (unsigned i = 0; i < testset.n; ++i) {
    printf("[%*d/%d]\t%-36s....  ", ndig, i+1, testset.n, testset.tests[i].name);
    fflush(stdout);
    if (run_test(&testset.tests[i])) {
      printf("Passed\n");
      ++npassed;
    }
  }
  bool pass = npassed == testset.n;
  if (pass) {
    printf("All %d tests passed for %s\n", testset.n, testset.name);
  } else {
    printf("\n%s: %d Passed | %d FAILED | %d Total\n",
           testset.name, npassed, testset.n - npassed, testset.n);
  }
  return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
