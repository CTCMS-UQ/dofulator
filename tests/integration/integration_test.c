#include <stdio.h>
#include <stdlib.h>

#include "integration_test.h"
#include "quaternion.h"
#include "fragment.h"

typedef struct RotationTest {
  const char* name;
  Quaternion q;
} RotationTest;

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
        printf("\t! Atom %zu expected %.16g DoF in %c, got %.16g\n",
               i, test->dof[i][d], (char[]){'x', 'y', 'z'}[d], dof[d]);
        result = false;
      }
      dof_atom += test->dof[i][d];
    }
    double dof_atom_actual = dofulator_get_dof_atom(ctx, i);
    if (!feql(dof_atom, dof_atom_actual)) {
      if (result) printf("FAILED\n");
      printf("\t! Atom %zu expected %.16g DoF total, got %.16g\n",
             i, dof_atom, dof_atom_actual);
      result = false;
    }
  }
  return result;
}

void rotate_system(Quaternion q, Test* test) {
  Quaternion q_conj = quat_conj(q);
  // Set up auxiliary system as ground truth.
  Dofulator ctx = dofulator_create(test->atoms.n);
  for (
    Bond* b = test->bonds.bonds;
    test->bonds.bonds && b < test->bonds.bonds + test->bonds.n;
    ++b
  ) {
    switch (test->test_type) {
      case RIGID:
        dofulator_build_rigid_fragment(ctx, *b);
        break;
      case FLEX:
        dofulator_add_rigid_bond(ctx, *b);
    }
  }
  dofulator_finalise_fragments(ctx);

  // For rigid tests, rotate the system before calculating to check
  // that rotated DoF calculation is correct
  if (test->test_type == RIGID) {
    for (size_t i = 0; i < test->atoms.n; ++i) {
      quat_rotate_vec(q, test->atoms.x[i]);
    }
  }
  dofulator_precalculate_rigid(ctx, test->atoms.mass, (const double(*)[3])test->atoms.x);
  dofulator_calculate(ctx, test->atoms.mass, (const double(*)[3])test->atoms.x);

  switch (test->test_type) {
    case RIGID:
      // Store result for comparison
      for (size_t i = 0; i < test->atoms.n; ++i) {
        dofulator_get_dof_atom_directional(ctx, i, test->dof[i]);
      }
      break;
    case FLEX:;
      // For flexible tests, rotate the DoF as for rigid fragments
      // for comparison against recalculating in new frame.
      double x[3] = {1., 0., 0.};
      double y[3] = {0., 1., 0.};
      double z[3] = {0., 0., 1.};
      quat_rotate_vec(q_conj, x);
      quat_rotate_vec(q_conj, y);
      quat_rotate_vec(q_conj, z);

      FragmentListIter frags = dofulator_get_semirigid_fragments(ctx);
      const Fragment* frag;
      while ((frag = fragmentlist_iter_next(&frags))) {
        for (uint32_t i = 0; i < frag->n_atoms; ++i) {
          AtomTag a = frag->atoms[i];
          // Sum up total DoF in each direction
          double d, rdx, rdy, rdz;
          test->dof[a][0] = test->dof[a][1] = test->dof[a][2] = 0.;
          for (size_t m = 0; m < frag->n_modes; ++m) {
            rdx = frag->dof[ 3*i    * 3*frag->n_atoms + m];
            rdy = frag->dof[(3*i+1) * 3*frag->n_atoms + m];
            rdz = frag->dof[(3*i+2) * 3*frag->n_atoms + m];
            d = rdx * x[0] + rdy * x[1] + rdz * x[2];
            test->dof[a][0] += d*d;
            d = rdx * y[0] + rdy * y[1] + rdz * y[2];
            test->dof[a][1] += d*d;
            d = rdx * z[0] + rdy * z[1] + rdz * z[2];
            test->dof[a][2] += d*d;
          }
        }
      }
      // Now apply rotation to system
      for (size_t i = 0; i < test->atoms.n; ++i) {
        quat_rotate_vec(q, test->atoms.x[i]);
      }
  }
}

bool run_test(Test* test) {
  bool result = true;
  Dofulator ctx = dofulator_create(test->atoms.n);
  if (!ctx) {
    printf("FAILED\n\tError creating dofulator context!\n");
    return false;
  }
  DofulatorResult err;
  bool init_success = false;
  const char* iter_msg = "";
  for (
    Bond* b = test->bonds.bonds;
    test->bonds.bonds && b < test->bonds.bonds + test->bonds.n;
    ++b
  ) {
    switch (test->test_type) {
      case RIGID:
        err = dofulator_build_rigid_fragment(ctx, *b);
        break;
      case FLEX:
        err = dofulator_add_rigid_bond(ctx, *b);
    }
    if (err) {
      result = false;
      goto cleanup;
    }
  }
  err = dofulator_finalise_fragments(ctx);
  if (err) {
    result = false;
    goto cleanup;
  }
  err = dofulator_precalculate_rigid(ctx, test->atoms.mass, (const double(*)[3])test->atoms.x);
  if (err) {
    result = false;
    goto cleanup;
  }
  init_success = true;

  iter_msg = "Initial calculation";
  err = dofulator_calculate(ctx, test->atoms.mass, (const double(*)[3])test->atoms.x);
  if (err) {
    result = false;
    goto cleanup;
  }
  result &= check_result(ctx, test, result);
  if (!result) goto cleanup;

  // Check for rotational invariance
  RotationTest rotations[] = {
    {
      .name = "Rotation invariance - 90 degrees around x axis",
      .q = {.x = sin(90.*DEG/2.), .y = 0, .z = 0, .w = cos(90*DEG/2.)},
    },
    {
      .name = "Rotation invariance - 2nd 90 degrees around x axis",
      .q = {.x = sin(90.*DEG/2.), .y = 0, .z = 0, .w = cos(90*DEG/2.)},
    },
    {
      .name = "Rotation invariance - 180 degrees around x axis",
      .q = {.x = sin(180.*DEG/2.), .y = 0, .z = 0, .w = cos(180*DEG/2.)},
    },
    {
      .name = "Rotation invariance - 30 degrees around z axis",
      .q = {.x = 0, .y = 0, .z = sin(30.*DEG/2.), .w = cos(30*DEG/2.)},
    },
    {
      .name = "Rotation invariance - Undo 30 degrees around z axis",
      .q = {.x = 0, .y = 0, .z = sin(-30.*DEG/2.), .w = cos(-30*DEG/2.)},
    },
    {
      .name = "Rotation invariance - 30 degrees around non-basis axis",
      .q = {.x = 0, .y = sin(30.*DEG/2.)*cos(30*DEG), .z = sin(30.*DEG/2.)*sin(30*DEG), .w = cos(30*DEG/2.)},
    },
    {
      .name = "Rotation invariance - Undo 30 degrees around non-basis axis",
      .q = {.x = 0, .y = sin(-30.*DEG/2.)*cos(30*DEG), .z = sin(-30.*DEG/2.)*sin(30*DEG), .w = cos(-30*DEG/2.)},
    },
  };

  for (RotationTest* rot = rotations; rot < rotations + sizeof(rotations)/sizeof(RotationTest); ++rot) {
    iter_msg = rot->name;
    rotate_system(rot->q, test);
    err = dofulator_calculate(ctx, test->atoms.mass, (const double(*)[3])test->atoms.x);
    if (err) {
      result = false;
      goto cleanup;
    }
    result &= check_result(ctx, test, result);
    if (!result) goto cleanup;
  }

cleanup:
  dofulator_destroy(&ctx);
  if (err) {
    printf("\t! Failure with error code: %d\n", err);
  }
  if (!init_success) {
    printf("\t! Failed to initialise context\n");
  } else if (!result) {
    printf("\t! Failure occurred on: %s\n", iter_msg);
  }
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
