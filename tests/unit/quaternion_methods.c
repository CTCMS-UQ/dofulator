#include "quaternion.h"
#include "unit_test.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {

  printf("Checking basic rotation \t\t\t... ");
  double a[3] = {1.0, 0.0, 0.0};
  double angle = 45*DEG;
  double b[3] = {cos(angle), sin(angle), 0.0};
  Quaternion q = {.x = 0, .y = 0, .z = sin(angle/2.), .w = cos(angle/2.)};
  quat_rotate_vec(q, a);
  bool result = true;
  result &= feql(a[0], b[0]);
  result &= feql(a[1], b[1]);
  result &= feql(a[2], b[2]);
  if (result) {
    printf("Passed\n");
  } else {
    printf("FAILED\n");
    return EXIT_FAILURE;
  }

  printf("Checking double rotation \t\t\t... ");
  a[0] = 1.; a[1] = a[2] = 0.;
  angle = 22.5*DEG;
  q = (Quaternion){.x = 0, .y = 0, .z = sin(angle/2.), .w = cos(angle/2.)};
  quat_rotate_vec(q, a);
  quat_rotate_vec(q, a);
  result = true;
  result &= feql(a[0], b[0]);
  result &= feql(a[1], b[1]);
  result &= feql(a[2], b[2]);
  if (result) {
    printf("Passed\n");
  } else {
    printf("FAILED\n");
    return EXIT_FAILURE;
  }

  printf("Checking 6x rotation \t\t\t\t... ");
  a[0] = 1.; a[1] = a[2] = 0.;
  b[0] = 1.; b[1] = b[2] = 0.;
  angle = 60*DEG;
  q = (Quaternion){.x = cos(45.*DEG)*sin(angle/2.), .y = sin(45.*DEG)*sin(angle/2.), .z = 0., .w = cos(angle/2.)};
  quat_rotate_vec(q, a);
  quat_rotate_vec(q, a);
  quat_rotate_vec(q, a);
  quat_rotate_vec(q, a);
  quat_rotate_vec(q, a);
  quat_rotate_vec(q, a);
  result = true;
  result &= feql(a[0], b[0]);
  result &= feql(a[1], b[1]);
  result &= feql(a[2], b[2]);
  if (result) {
    printf("Passed\n");
  } else {
    printf("FAILED\n");
    return EXIT_FAILURE;
  }

  printf("Checking quat_from_closest_arc, simple case\t... ");
  a[0] = a[1] = a[2] = 1.0;
  b[0] = b[1] = 1.0 / sqrt(2.);
  b[2] = 0.;
  q = quat_from_closest_arc(a, b);
  quat_rotate_vec(q, a);
  result = true;
  result &= feql(a[2], 0.);
  result &= feql(a[0], a[1]);
  result &= feql(vec_dot(a, a), 3.);
  if (result) {
    printf("Passed\n");
  } else {
    printf("FAILED\n");
    return EXIT_FAILURE;
  }

  printf("Checking quat_from_closest_arc, 180 deg\t\t... ");
  a[0] = -1.;
  a[1] = -1.;
  a[2] = 0.;
  q = quat_from_closest_arc(a, b);
  quat_rotate_vec(q, a);
  result = true;
  result &= feql(a[2], 0.);
  result &= feql(a[0], a[1]);
  result &= feql(vec_dot(a, a), 2.);
  result &= feql(a[0] / sqrt(2.), b[0]);
  result &= feql(a[1] / sqrt(2.), b[1]);
  if (result) {
    printf("Passed\n");
  } else {
    printf("FAILED\n");
    return EXIT_FAILURE;
  }

  printf("Checking compound rotation\t\t\t... ");
  a[0] = 0.;
  a[1] = 1.;
  a[2] = 0.;

  b[0] = 0.;
  b[1] = 0.;
  b[2] = 1.;

  // Want to rotate the frame defined by {aa, bb} onto {a, b}
  double aa[3] = {1., 0., 0.};
  double bb[3] = {0., 1.0/sqrt(2.), -1./sqrt(2.)};

  // Rotation for aa onto a
  Quaternion qa = quat_from_closest_arc(aa, a);
  quat_rotate_vec(qa, bb);

  // Rotation around a to get bb onto b
  Quaternion qb = quat_from_vec(a);
  double cos_angle = vec_dot(bb, b);
  if (fabs(1. + cos_angle) > 100.*DBL_EPSILON) {
    double sin_angle = sqrt(1 - cos_angle*cos_angle);
    qb.x *= sin_angle;
    qb.y *= sin_angle;
    qb.z *= sin_angle;
  }
  qb.w = 1. + cos_angle;
  qb = quat_normalize(qb);

  q = quat_mul(qb, qa);

  // Check that the compound rotation works
  bb[0] = 0.;
  bb[1] = 1.0/sqrt(2.);
  bb[2] = -1.0/sqrt(2.);
  quat_rotate_vec(q, aa);
  quat_rotate_vec(q, bb);
  result = true;
  result &= feql(aa[0], a[0]);
  result &= feql(aa[1], a[1]);
  result &= feql(aa[2], a[2]);
  result &= feql(bb[0], b[0]);
  result &= feql(bb[1], b[1]);
  result &= feql(bb[2], b[2]);
  if (result) {
    printf("Passed\n");
  } else {
    printf("FAILED\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
