#include "quaternion.h"
#include "unit_test.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {

  printf("Checking quat_from_closest_arc, simple case\t... ");
  double a[3] = {1.0, 1.0, 1.0};
  double b[3] = {1.0 / sqrt(2.), 1.0 / sqrt(2.), 0.0};
  Quaternion q = quat_from_closest_arc(a, b);
  quat_rotate_vec(q, a);
  assert(feql(a[2], 0.));
  assert(feql(a[0], a[1]));
  assert(feql(vec_dot(a, a), 3.));
  printf("Passed\n");

  printf("Checking quat_from_closest_arc, 180 deg\t\t... ");
  a[0] = -1.;
  a[1] = -1.;
  a[2] = 0.;
  q = quat_from_closest_arc(a, b);
  quat_rotate_vec(q, a);
  assert(feql(a[2], 0.));
  assert(feql(a[0], a[1]));
  assert(feql(vec_dot(a, a), 2.));
  assert(feql(a[0] / sqrt(2.), b[0]));
  assert(feql(a[1] / sqrt(2.), b[1]));
  printf("Passed\n");

  printf("Checking compound rotation\t\t\t... ");
  a[0] = 0.;
  a[1] = 1.;
  a[2] = 0.;

  b[0] = 0.;
  b[1] = 0.;
  b[2] = 1.;

  double aa[3] = {1., 0., 0.};
  double bb[3] = {0., 1.0/sqrt(2.), -1./sqrt(2.)};

  Quaternion qa = quat_from_closest_arc(aa, a);
  quat_rotate_vec(qa, bb);

  Quaternion qb = quat_from_vec(a);
  double cos_angle = vec_dot(bb, b);
  double sin_angle = sqrt(1 - cos_angle*cos_angle);
  qb.x *= sin_angle;
  qb.y *= sin_angle;
  qb.z *= sin_angle;
  qb.w = sqrt(vec_dot(bb, bb)) + cos_angle;
  qb = quat_normalize(qb);

  q = quat_mul(qb, qa);

  bb[0] = 0.;
  bb[1] = 1.0/sqrt(2.);
  bb[2] = -1.0/sqrt(2.);
  quat_rotate_vec(q, aa);
  quat_rotate_vec(q, bb);
  assert(feql(aa[0], a[0]));
  assert(feql(aa[1], a[1]));
  assert(feql(aa[2], a[2]));
  assert(feql(bb[0], b[0]));
  assert(feql(bb[1], b[1]));
  assert(feql(bb[2], b[2]));
  printf("Passed\n");
  
  return EXIT_SUCCESS;
}
