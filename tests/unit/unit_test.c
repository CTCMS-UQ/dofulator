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
  Fragment frag;
  switch (test->test_type) {
    case RIGID:
      frag = fragment_create_rigid(test->atoms);
      break;
    case FLEX:
      frag = fragment_create_semirigid(test->atoms, test->bonds);
      break;
    default:
      return false;
  }
  if (!frag) {
    printf("FAILED\n\tError creating rigid body!\n");
    return false;
  }

  double dof_total = 0.0;
  for (unsigned i = 0; i < test->atoms.n; ++i) {
    double dof_atom = 0.0;
    for (unsigned d = 0; d < 3; ++d) {
      double dir[] = {0.0, 0.0, 0.0};
      dir[d] = 1.0;
      double dof = fragment_dof_atom_dir(frag, i, dir);
      if (!feql(dof, test->dof[i][d])) {
        if (result) printf("FAILED\n");
        printf("\t! Atom %d expected %.16g DoF in %c, got %.16g\n",
               i, test->dof[i][d], (char[]){'x', 'y', 'z'}[d], dof);
        result = false;
      }
      dof_atom += test->dof[i][d];
      dof_total += test->dof[i][d];
    }
    double dof = fragment_dof_atom(frag, i);
    if (!feql(dof_atom, dof)) {
      if (result) printf("FAILED\n");
      printf("\t! Atom %d expected %.16g DoF, got %.16g\n",
             i, dof_atom, dof);
      result = false;
    }
  }
  double dof = fragment_dof(frag);
  if (!feql(dof_total, dof)) {
    if (result) printf("FAILED\n");
    printf("\t! Expected %.16g DoF total, got %.16g\n",
           dof_total, dof);
    result = false;
  }

  fragment_destroy(frag);
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
