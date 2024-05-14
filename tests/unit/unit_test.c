#include <stdio.h>
#include <stdlib.h>

#include "unit_test.h"
#include "quaternion.h"

static unsigned digits(unsigned n) {
  unsigned d = 1;
  while (n >= 10) {
    n /= 10;
    ++d;
  }
  return d;
}

bool check_result(Dofulator ctx, Test* test, bool current_result) {
  bool result = current_result;
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
    }
    double dof_atom_actual = dofulator_get_dof_atom(ctx, i);
    if (!feql(dof_atom, dof_atom_actual)) {
      if (result) printf("FAILED\n");
      printf("\t! Atom %ld expected %.16g DoF total, got %.16g\n",
             i, dof_atom, dof_atom_actual);
      result = false;
    }
  }
  return result;
}

void rotate_system(Quaternion q, Test* test) {
  double x[3] = {1., 0., 0.};
  double y[3] = {0., 1., 0.};
  double z[3] = {0., 0., 1.};
  double dof[3];
  quat_rotate_vec(q, x);
  quat_rotate_vec(q, y);
  quat_rotate_vec(q, z);
  for (size_t i = 0; i < test->atoms.n; ++i) {
    quat_rotate_vec(q, test->atoms.x[i]);
    dof[0] = test->dof[i][0]*x[0]*x[0] + test->dof[i][1]*y[0]*y[0] + test->dof[i][2]*z[0]*z[0];
    dof[1] = test->dof[i][0]*x[1]*x[1] + test->dof[i][1]*y[1]*y[1] + test->dof[i][2]*z[1]*z[1];
    dof[2] = test->dof[i][0]*x[2]*x[2] + test->dof[i][1]*y[2]*y[2] + test->dof[i][2]*z[2]*z[2];
    test->dof[i][0] = dof[0];
    test->dof[i][1] = dof[1];
    test->dof[i][2] = dof[2];
  }
}

bool run_test(Test* test) {
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
  result &= check_result(ctx, test, result);

  // Check for rotational invariance
  double angle = 90. * DEG;
  Quaternion q = {.x = sin(angle/2.), .y = 0, .z = 0, .w = cos(angle/2.)};
  rotate_system(q, test);
  dofulator_calculate(ctx, test->atoms.mass, test->atoms.x);
  result &= check_result(ctx, test, result);

  // Cover 180 degree rotation
  rotate_system(q, test);
  dofulator_calculate(ctx, test->atoms.mass, test->atoms.x);
  result &= check_result(ctx, test, result);

  // Non-90 degree rotation
  angle = 30. * DEG;
  q = (Quaternion){.x = 0, .y = cos(30*DEG)*sin(angle/2.), .z = sin(30*DEG)*sin(angle/2.), .w = cos(angle/2.)};
  rotate_system(q, test);
  dofulator_calculate(ctx, test->atoms.mass, test->atoms.x);
  result &= check_result(ctx, test, result);

  dofulator_destroy(&ctx);
  return result;
}

int run_testset(TestSet testset) {
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
