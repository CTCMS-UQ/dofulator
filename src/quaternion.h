#ifndef DOFULATOR_QUAT_H
#define DOFULATOR_QUAT_H

#include <assert.h>
#include <float.h>
#include <math.h>

#include "compat.h"
#include "vec3.h"

typedef struct Quaternion {
  double x, y, z, w;
} Quaternion;

// Convert vector to quaternion (w = 0)
static inline Quaternion quat_from_vec(const double vec[3]) {
  return (Quaternion){ .x = vec[0], .y = vec[1], .z = vec[2], .w = 0. };
}

// Normalize a quaternion
static inline Quaternion quat_normalize(Quaternion q) {
  double len = sqrt(q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w);
  q.x /= len;
  q.y /= len;
  q.z /= len;
  q.w /= len;
  return q;
}

// Multiplication of two quaternions (not commutative)
static inline Quaternion quat_mul(Quaternion q1, Quaternion q2) {
  return (Quaternion){
    .x = q1.x*q2.w + q1.w*q2.x + q1.y*q2.z - q1.z*q2.y,
    .y = q1.y*q2.w + q1.w*q2.y + q1.z*q2.x - q1.x*q2.z,
    .z = q1.z*q2.w + q1.w*q2.z + q1.x*q2.y - q1.y*q2.x,
    .w = q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z,
  };
}

// Conjugate of a quaternion
static inline Quaternion quat_conj(Quaternion q) {
  q.x *= -1;
  q.y *= -1;
  q.z *= -1;
  return q;
}

// Apply the rotation of q to v in-place
static inline void quat_rotate_vec(Quaternion q, double v[3]) {
  q = quat_mul(quat_mul(q, quat_from_vec(v)), quat_conj(q));
  v[0] = q.x; v[1] = q.y; v[2] = q.z;
}

// Identity quaternion
static inline Quaternion quat_identity(void) {
  return (Quaternion){.x = 0., .y = 0., .z = 0., .w = 1.};
}

// Find quaternion which rotates v1 onto v2, assuming v2 is normalized
static inline Quaternion quat_from_closest_arc(const double NOALIAS_ARR(v1, 3), const double NOALIAS_ARR(v2, 3)) {
  assert(vec_dot(v2, v2) - 1. < 100. * DBL_EPSILON);

  // For non-normalised vectors, general formula is:
  // q.{x, y, z} = (v1 x v2)
  // q.w = |v1||v2| + v1 . v2
  // see http://www.euclideanspace.com/maths/algebra/vectors/angleBetween/minorlogic.htm
  // Adapted to handle 180 degree angles

  double w = sqrt(vec_dot(v1, v1)) + vec_dot(v1, v2);
  double c[3];
  if (w > 100. * DBL_EPSILON) {
    vec_cross(v1, v2, c);
  } else {
    // 180 degree angle - find vector perpendicular to v2 as axis for 180deg rotation
    // see https://math.stackexchange.com/a/4112622
    w = 0.;
    c[0] = copysign(v2[2], v2[0]);
    c[1] = copysign(v2[2], v2[1]);
    c[2] = copysign(fabs(v2[0]) + fabs(v2[1]), v2[2]);
  }
  return quat_normalize((Quaternion){
    .x = c[0], .y = c[1], .z = c[2], .w = w
  });
}

#endif
