#ifndef DOFULATOR_UNIT_TEST_H
#define DOFULATOR_UNIT_TEST_H

#include <float.h>
#include <math.h>
#include <stdbool.h>

#include "dofulator.h"

#define PI 3.14159265378979
#define DEG (PI/180.0)

#define DBLS (double*)(double[])
#define ntests(testset) (sizeof(testset)/sizeof(*testset))

inline static bool feql(double a, double b) {
  return fabs(a - b) < 100*DBL_EPSILON;
}

typedef struct Test {
  enum {RIGID, FLEX} test_type;
  char* name;
  AtomList atoms;
  BondList bonds;
  double (*dof)[3]; // Expected DoF in each direction for each atom
} Test;

typedef struct TestSet {
  char* name;
  unsigned n;
  Test* tests;
} TestSet;

bool run_test(const Test* test);
int run_testset(const TestSet testset);

#endif
