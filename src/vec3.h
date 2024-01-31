#ifndef VEC3_H
#define VEC3_H

typedef struct Vec3 {
  double x;
  double y;
  double z;
} Vec3;

// Returns l - r
inline static Vec3 vec_sub(double l[3], double r[3]) {
  return (Vec3){
    .x = l[0] - r[0],
    .y = l[1] - r[1],
    .z = l[2] - r[2]
  };
}

#endif
