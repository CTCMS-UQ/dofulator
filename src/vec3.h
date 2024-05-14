#ifndef DOFULATOR_VEC3_H
#define DOFULATOR_VEC3_H

#include <math.h>

inline static void vec_sub(const double l[restrict 3], const double r[restrict 3], double out[restrict 3]) {
  out[0] = l[0] - r[0];
  out[1] = l[1] - r[1];
  out[2] = l[2] - r[2];
}

static inline void vec_cross(const double v1[restrict 3], const double v2[restrict 3], double c[restrict 3]) {
  c[0] = v1[1] * v2[2] - v1[2] * v2[1];
  c[1] = v1[2] * v2[0] - v1[0] * v2[2];
  c[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

static inline double vec_dot(const double v1[3], const double v2[3]) {
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

static inline void vec_normalize(double v[3]) {
  double len = sqrt(vec_dot(v, v));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

#endif
