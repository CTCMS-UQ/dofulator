#ifndef FRAGMENT_H
#define FRAGMENT_H

#include "dofulator.h"

typedef struct Fragment_t {
  unsigned natoms;  // Number of atoms in the fragment
  unsigned nmodes;  // Maximum number of modes. Some could have zero inertia, in
                    // which case they are not a real DoF.

  double *mC;       // M^1/2 C, where C = constraint matrix, M = diagonal per-atom,
                    // per-direction mass matrix. C^T C = lab-frame inertia of
                    // internal modes (Itot) [3*natoms][nmodes]

  double* mCR;      // Stores M^1/2 C R, where R = eigenvectors of C^T M C.
                    // R^T Itot R == mCR^T mCR is diagonal. [3*natoms][nmodes]

  double* Imodal;   // Total modal inertia == diag(mCR^T mCR), summed over sets of 3
                    // for 3 cartesian directions of each atom [natoms]

  double R[];       // Eigenvectors of Itot == C^T M C == mC^T mC. [nmodes][nmodes]

} Fragment_t;

Fragment_t* fragment_eval(Fragment_t* frag);

static inline Fragment_t* fragment_alloc(AtomList atoms, unsigned max_modes);

#endif
