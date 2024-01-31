#ifndef FRAGMENT_H
#define FRAGMENT_H

#include "dofulator.h"

typedef struct Fragment_t {
  AtomList atoms;   // Atoms part of the rigid/semi-rigid fragment. Not owned, must outlive the fragment.
  unsigned n_modes; // Number of modes (upper bound)
  double* M;        // per-atom mass matrix
  double* C;        // Constraint matrix
  double* R;        // Eigenvectors of Itot. R^T R == 1, R^T Itot R == diagonal
  double* Iatom;    // Modal inertia of each atom == R^T C^T M_i C R (length 3*atoms.n * 3*atoms.n)
  double* Imodal;   // Total modal inertia == R^T C^T M C R (length atoms.n)
  double buf[];     // Data storage block
} Fragment_t;

Fragment_t* fragment_eval(Fragment_t* frag);

static inline Fragment_t* fragment_alloc(AtomList atoms);

#endif
