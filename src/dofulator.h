#ifndef LIBDOFULATOR_H
#define LIBDOFULATOR_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef AtomTag
#define AtomTag size_t
#endif

// TODO: This should be left at 1 for now since batch operations
//       aren't implemented properly yet.
#ifndef DOFULATOR_DEFAULT_BATCH_SIZE 
#define DOFULATOR_DEFAULT_BATCH_SIZE 1
#endif

typedef struct Dofulator* Dofulator;
typedef struct Fragment Fragment;

typedef struct Bond {
  AtomTag ai;
  AtomTag aj;
} Bond;

typedef struct FragmentListIter {
  size_t n_fragments;         // Number of fragments to iterate over
  size_t idx;                 // Current index
  const Fragment* fragments;  // Read-only view of fragments. Invalidated if FragmentList changes.
} FragmentListIter;


typedef struct AtomListView {
  size_t n_atoms;         // Number of atoms.
  const AtomTag* atoms;   // Atom indices.
                          // Read-only, but may be changed/invalidated if
                          // fragments are rebuilt. Should be copied out
                          // if persistence is needed.
} AtomListView;


typedef struct PBC {
  enum { PBC_NONE = 0, PBC_TRI, PBC_ORTHO} typ;
  union {
    struct {
      double lx, _pad_yzx[3], ly, _pad_zxy[3], lz;
    };
    struct {
      double ax, _pad_a[2];
      double bx, by, _pad_b;
      double cx, cy, cz;
    };
    struct {
      double a[3], b[3], c[3];
    };
  };
} PBC;


Dofulator dofulator_create(AtomTag n_atoms);
void dofulator_destroy(Dofulator* ctx);
void dofulator_set_pbc(Dofulator ctx, PBC pbc);

void dofulator_add_rigid_bond(Dofulator ctx, Bond b);
void dofulator_build_rigid_fragment(Dofulator ctx, Bond b);
void dofulator_finalise_fragments(Dofulator ctx);

void dofulator_precalculate_rigid(Dofulator ctx, const double* mass, const double x[][3]);
void dofulator_calculate(Dofulator ctx, const double* mass, const double x[][3]);

void dofulator_get_dof_atom_directional(const struct Dofulator* ctx, AtomTag atom_idx, double dof[3]);
double dofulator_get_dof_atom(const struct Dofulator* ctx, AtomTag atom_idx);

FragmentListIter dofulator_get_rigid_fragments(const struct Dofulator* ctx);
FragmentListIter dofulator_get_semirigid_fragments(const struct Dofulator* ctx);

const Fragment* fragmentlist_iter_next(FragmentListIter* iter);

AtomListView fragment_get_atoms(const Fragment* frag);

#ifdef __cplusplus
}
#endif

#endif
