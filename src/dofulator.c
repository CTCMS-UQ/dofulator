#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef NDEBUG
#include <stdio.h>
#endif

#include "cblas_lapacke.h"
#include "dofulator.h"
#include "fragment.h"
#include "quaternion.h"
#include "vec3.h"

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

/*******************************************************************************
 * Main context object
*/
struct Dofulator {
  AtomTag n_atoms;                // Number of atoms in frag_map
  FragIndex* frag_map;            // Index of fragment corresponding to given atom
  AtomTag* predecessors;          // Predecessor for each atom (self for fragment root nodes)
  size_t* atom_frag_idx;          // Index of the atom within the fragment. Gives corresponding rows in Jacobian
  size_t* n_semirigid;            // Number of semirigid fragments indexed by (number of atoms - 1)
  size_t* n_rigid;                // Number of rigid fragments indexed by number of atoms
  size_t max_semirigid_atoms;     // Maximum number of atoms in any semirigid fragment
  size_t max_rigid_atoms;         // Maximum number of atoms in any rigid fragment
  size_t max_K_size;              // Maximum storage space needed by a K matrix
  size_t batch_size;              // Minimum number of fragments to process at once. TODO: overhaul working memory allocation to enable proper batching.
  FragmentList semirigid_frags;   // Fragment list, partitioned by Jacobian size (number of atoms)
  FragmentList rigid_frags;       // Fragment list, partitioned by Jacobian size (number of atoms)
  RefFrame* rigid_ref_frames;     // Reference frame for each rigid fragment
  double** dof_buf_semirigid;     // Buffer for each fragment to store its dof result, partitioned by number of atoms.
  double** dof_total_buf_semirigid;  // Buffer for each fragment to store its dof_total result, partitioned by numberof atoms.
  double** dof_buf_rigid;         // Buffer for each fragment to store its dof result, partitioned by number of atoms.
  double** dof_total_buf_rigid;   // Buffer for each fragment to store its dof_total result, partitioned by numberof atoms.
  size_t working_buf_size;        // Size of working memory
  double* working_buf;            // Working memory for building Jacobians
  double* svd_working;            // Working memory for SVD/eigenvector calculation (reallocated as necessary)
  double null_space_thresh;       // Fraction of the maximum singular value below which singular values are treated as 0.
                                  // 0.0 = use DBL_EPSILON * max dimension of loop closure matrix.
  lapack_int svd_working_size;    // Size of SVD/eigenvector working memory
  PBC pbc;
};


/*******************************************************************************
 * Create a dofulator context with space for `n_atoms` atoms.
 * Returns NULL on allocation failure.
*/
Dofulator dofulator_create(AtomTag n_atoms) {
  struct Dofulator ctx = {
    .n_atoms = n_atoms,
    .batch_size = DOFULATOR_DEFAULT_BATCH_SIZE,
    .frag_map = (FragIndex*)calloc(n_atoms, sizeof(FragIndex)),
    .predecessors = (AtomTag*)malloc(sizeof(AtomTag) * n_atoms),
    .atom_frag_idx = (size_t*)malloc(sizeof(size_t) * n_atoms),
    .null_space_thresh = 0.0,   // Unlikely to need changes, but expose just in case
  };
  for (AtomTag i = 0; i < n_atoms; ++i) {
    ctx.predecessors[i] = i;
  }
  Dofulator out = malloc(sizeof(ctx));
  *out = ctx;
  return out;
}


/*******************************************************************************
 * Clean up a dofulator context.
*/
void dofulator_destroy(Dofulator* ctx) {
  Dofulator self = *ctx;
  *ctx = NULL;
  if (!self) return;
  if (self->frag_map) free(self->frag_map);
  if (self->predecessors) free(self->predecessors);
  if (self->atom_frag_idx) free(self->atom_frag_idx);
  if (self->n_semirigid) free(self->n_semirigid);
  if (self->n_rigid) free(self->n_rigid);

  if (self->dof_buf_semirigid) {
    for (double** b = self->dof_buf_semirigid; b < self->dof_buf_semirigid + self->max_semirigid_atoms; ++b) {
      if (*b) free(*b);
    }
    free(self->dof_buf_semirigid);
  }
  if (self->dof_total_buf_semirigid) {
    for (double** b = self->dof_total_buf_semirigid; b < self->dof_total_buf_semirigid + self->max_semirigid_atoms; ++b) {
      if (*b) free(*b);
    }
    free(self->dof_total_buf_semirigid);
  }
  self->max_semirigid_atoms = 0;
  self->max_K_size = 0;

  if (self->dof_buf_rigid) {
    for (double** b = self->dof_buf_rigid; b < self->dof_buf_rigid + self->max_rigid_atoms; ++b) {
      if (*b) free(*b);
    }
    free(self->dof_buf_rigid);
  }
  if (self->dof_total_buf_rigid) {
    for (double** b = self->dof_total_buf_rigid; b < self->dof_total_buf_rigid + self->max_rigid_atoms; ++b) {
      if (*b) free(*b);
    }
    free(self->dof_total_buf_rigid);
  }
  self->max_rigid_atoms = 0;

  if (self->working_buf) free(self->working_buf);
  if (self->svd_working) free(self->svd_working);
  fragmentlist_destroy(&self->semirigid_frags);
  fragmentlist_destroy(&self->rigid_frags);

  self->n_atoms = 0;
  self->working_buf_size = 0;
  self->svd_working_size = 0;
}

/*******************************************************************************
 * Setter for pbc via opaque handle
*/
void dofulator_set_pbc(Dofulator ctx, PBC pbc) {
  ctx->pbc = pbc;
}

/*******************************************************************************
 * Setter for null_space_thresh via opaque handle
*/
void dofulator_set_null_space_thresh(Dofulator ctx, double thresh) {
  thresh = fabs(thresh);
  ctx->null_space_thresh = thresh < 1. ? thresh : 1.;
}

/*******************************************************************************
 * Getter for null_space_thresh via opaque handle
*/
double dofulator_get_null_space_thresh(const struct Dofulator* ctx) {
  return ctx->null_space_thresh;
}


/*******************************************************************************
 * Remap vector `r` to be within +/- 0.5*[lx, ly, lz] if necessary
*/
static inline void pbc_wrap_shortest(const PBC* pbc, double r[3]) {
  switch (pbc->typ) {
    case PBC_NONE:
      return;

    case PBC_TRI:
      // z
      if (pbc->lz != 0.) {
        if (r[2] > pbc->lz) {
          int n = (int)(r[2]/pbc->lz);
          r[0] -= n * pbc->cx;
          r[1] -= n * pbc->cy;
          r[2] -= n * pbc->cz;
        }
        while (r[2] > pbc->lz/2.) {
          r[0] -= pbc->cx;
          r[1] -= pbc->cy;
          r[2] -= pbc->cz;
        }
        if (-r[2] > pbc->lz) {
          int n = (int)(-r[2]/pbc->lz);
          r[0] += n * pbc->cx;
          r[1] += n * pbc->cy;
          r[2] += n * pbc->cz;
        }
        while (r[2] < -pbc->lz/2.) {
          r[0] += pbc->cx;
          r[1] += pbc->cy;
          r[2] += pbc->cz;
        }
      }

      // y
      if (pbc->ly != 0.) {
        if (r[1] > pbc->ly) {
          int n = (int)(r[1]/pbc->ly);
          r[0] -= n * pbc->bx;
          r[1] -= n * pbc->by;
        }
        while (r[1] > pbc->ly/2.) {
          r[0] -= pbc->bx;
          r[1] -= pbc->by;
        }
        if (-r[1] > pbc->ly) {
          int n = (int)(-r[1]/pbc->ly);
          r[0] += n * pbc->bx;
          r[1] += n * pbc->by;
        }
        while (r[1] < -pbc->ly/2.) {
          r[0] += pbc->bx;
          r[1] += pbc->by;
        }
      }

      // x
      if (pbc->lz != 0.) {
        if (r[0] > pbc->lx) {
          int n = (int)(r[0]/pbc->lx);
          r[0] -= n * pbc->ax;
        }
        while (r[0] > pbc->lx/2.) {
          r[0] -= pbc->ax;
        }
        if (-r[0] > pbc->lz) {
          int n = (int)(-r[0]/pbc->lx);
          r[0] += n * pbc->ax;
        }
        while (r[0] < -pbc->lx/2.) {
          r[0] += pbc->ax;
        }
      }
      break;

    case PBC_ORTHO:
      // x
      if (pbc->lx != 0.) {
        if (r[0] > pbc->lx) {
          r[0] -= ((int)(r[0]/pbc->lx)) * pbc->lx;
        }
        while (r[0] >  pbc->lx/2.) r[0] -= pbc->lx;
        if (-r[0] > pbc->lx) {
          r[0] += ((int)(-r[0]/pbc->lx)) * pbc->lx;
        }
        while (r[0] < -pbc->lx/2.) r[0] += pbc->lx;
      }

      // y
      if (pbc->ly != 0.) {
        if (r[1] > pbc->ly) {
          r[1] -= ((int)(r[1]/pbc->ly)) * pbc->ly;
        }
        while (r[1] >  pbc->ly/2.) r[1] -= pbc->ly;
        if (-r[1] > pbc->ly) {
          r[1] += ((int)(-r[1]/pbc->ly)) * pbc->ly;
        }
        while (r[1] < -pbc->ly/2.) r[1] += pbc->ly;
      }

      // z
      if (pbc->lz != 0.) {
        if (r[2] > pbc->lz) {
          r[2] -= ((int)(r[2]/pbc->lz)) * pbc->lz;
        }
        while (r[2] >  pbc->lz/2.) r[2] -= pbc->lz;
        if (-r[2] > pbc->lz) {
          r[2] += ((int)(-r[2]/pbc->lz)) * pbc->lz;
        }
        while (r[2] < -pbc->lz/2.) r[2] += pbc->lz;
      }
      break;
  }
}


/*******************************************************************************
 * Set the batch size and reallocate working memory
 *
 * Working memory:
 * Max. Jacobian size = MAX((max_semirigid_atoms * 3)^2, (max_rigid_atoms*3 * 6))
 * Max. Eigenvector size = MAX((max_semirigid_atoms * 3)^2, 6^2)
 *
 * Solve process:
 * Build J (in mJQ for loops, or in mJ for trees/rigid)
 * (loops only)
 *    build K from J (can use mJ scratch space for result, VT scratch for J_j - J_i, Q scratch for T^T) -> batched dgemm for each T^T (J_j - J_i)
 *    find null(K) (store in VT scratch space, singular values in dof_total)
 *      -> dgesvd in loop over batch. Pass NULL for U.
 *      -> Also needs WORK storage (get size from dgesvd call, can allocate for single fragment and re-use)
 *    project J onto null(K) (store in mJ scratch space) -> batched dgemm
 * Calculate sqrt(m) J (directly in scratch mJ space) -> dscal (batched?) modifies in-place
 * Calculate mJ^T mJ (store in Q scratch space) -> batched dgemm
 * Calculate eigenvectors (result replaces Q) -> dsyev in loop over batch, eigenvalues in Itotal
 * Calculate mJQ = mJ (mJ scratch space) * Q (Q scratch space) -> batched dgemm
 * Calculate dof[i][m] = mJQ[i][m] mJQ[i][m] / Itotal[m]
 *
 * Storage requirements per batch:
 * - mJ space (Jacobian size, or [n_loops * 3*n_atoms] for loop closures)
 * - Q  space (Eigenvector size)
 * - VT space (Jacobian size - semi-rigid only)
 * - WORK space for SVD (pick largest required by any loop-closing fragments)
*/
static inline DofulatorResult dofulator_update_batch_size(Dofulator ctx, size_t batch_size) {
  ctx->batch_size = batch_size;
  size_t semirigid_size = 3 * ctx->max_semirigid_atoms * 3 * ctx->max_semirigid_atoms;
  size_t max_jacobian = MAX(MAX( 3 * ctx->max_rigid_atoms * 6, semirigid_size), ctx->max_K_size);
  size_t max_eigenvector = MAX(36, semirigid_size);
  size_t max_scratch_size_single = max_jacobian + max_eigenvector + semirigid_size;
  ctx->working_buf_size = batch_size * max_scratch_size_single;
  double* new_ptr = realloc(ctx->working_buf, sizeof(double) * ctx->working_buf_size);
  if (unlikely(!new_ptr)) {
    return DOF_ALLOC_FAILURE;
  }
  ctx->working_buf = new_ptr;
  return DOF_SUCCESS;
}


/*******************************************************************************
 * Get the index of the semi-rigid fragment which contains the atom with index `atom_idx`.
 * May mutate `ctx->frag_map` if it points to an invalidated (merged) fragment.
*/
static inline size_t dofulator_get_semirigid_idx(Dofulator ctx, AtomTag atom_idx) {
  return fragmentlist_get_fragment_idx(ctx->semirigid_frags, ctx->frag_map, atom_idx);
}

/*******************************************************************************
 * Get the index of the rigid fragment which contains the atom with index `atom_idx`.
 * May mutate `ctx->frag_map` if it points to an invalidated (merged) fragment.
*/
static inline size_t dofulator_get_rigid_idx(Dofulator ctx, AtomTag atom_idx) {
  return fragmentlist_get_fragment_idx(ctx->rigid_frags, ctx->frag_map, atom_idx);
}


/*******************************************************************************
 * Update list of semirigid fragments to account for a rigid bond.
*/
DofulatorResult dofulator_add_rigid_bond(Dofulator ctx, Bond b) {
  if (unlikely(b.ai < 0 || b.ai >= ctx->n_atoms || b.aj < 0 || b.aj >= ctx->n_atoms)) {
    return DOF_INDEX_ERROR;
  }
  if (unlikely(!ctx)) { return DOF_UNINITIALISED; }
  if (unlikely(b.ai == b.aj)) { return DOF_SUCCESS; }

  if (ctx->frag_map[b.ai].has_frag && ctx->frag_map[b.aj].has_frag) {
    if (unlikely(ctx->frag_map[b.ai].rigid || ctx->frag_map[b.aj].rigid)) {
      return DOF_MIXED_RIGID_SEMIRIGID;
    }
    size_t idx_i = dofulator_get_semirigid_idx(ctx, b.ai);
    size_t idx_j = dofulator_get_semirigid_idx(ctx, b.aj);
    // Both already in a fragment
    if (idx_i == idx_j) {
      // Already in same fragment, so check if duplicate bond
      if (ctx->predecessors[b.ai] == b.aj || ctx->predecessors[b.aj] == b.ai) {
        return DOF_SUCCESS;
      }
      // Not duplicate, so this is a loop closure
      Fragment* frag = &ctx->semirigid_frags.fragments[idx_i];
      Bond* new_ptr = realloc(frag->loop_closures, sizeof(Bond) * (frag->n_loops + 1));
      if (unlikely(!new_ptr)) { return DOF_ALLOC_FAILURE; }
      frag->loop_closures = new_ptr;
      frag->loop_closures[frag->n_loops++] = b;
      return DOF_SUCCESS; // Return here so predecessors aren't changed
    } else {
      // Two different fragments, so join them
      Fragment* frag_i = &ctx->semirigid_frags.fragments[idx_i];
      Fragment* frag_j = &ctx->semirigid_frags.fragments[idx_j];
      assert(frag_i->root_atom != frag_j->root_atom);
      frag_i->n_atoms += frag_j->n_atoms;
      ctx->frag_map[b.aj].idx = idx_i;
      // Invalidate frag_j and point it to frag_i
      // Any other atoms in it will be updated later through this link
      // via `dofulator_get_semirigid_idx(ctx, idx_j)`
      frag_j->invalid = true;
      frag_j->root_atom = (AtomTag)idx_i;
      if (frag_j->n_loops > 0) {
        // Transfer loop closures across
        Bond* new_ptr = realloc(frag_i->loop_closures, sizeof(Bond) * (frag_i->n_loops + frag_j->n_loops));
        if (unlikely(!new_ptr)) { return DOF_ALLOC_FAILURE; }
        frag_i->loop_closures = new_ptr;
        memcpy(&frag_i->loop_closures[frag_i->n_loops], frag_j->loop_closures, sizeof(Bond) * frag_j->n_loops);
        frag_i->n_loops += frag_j->n_loops;
        frag_j->n_loops = 0;
        free(frag_j->loop_closures);
        frag_j->loop_closures = NULL;
      }
      // Reorient tree of frag_j to correct predecessors
      AtomTag new_pred = b.ai;
      AtomTag atom = b.aj;
      AtomTag old_pred = ctx->predecessors[atom];
       while (atom != new_pred) {
        ctx->predecessors[atom] = new_pred;
        new_pred = atom;
        atom = old_pred;
        old_pred = ctx->predecessors[atom];
      }
    }

  } else if (ctx->frag_map[b.ai].has_frag) {
    if (unlikely(ctx->frag_map[b.ai].rigid)) { return DOF_MIXED_RIGID_SEMIRIGID; }
    // Add j to i's fragment
    size_t idx_i = dofulator_get_semirigid_idx(ctx, b.ai);
    Fragment* frag = &ctx->semirigid_frags.fragments[idx_i];
    frag->n_atoms++;
    ctx->frag_map[b.aj] = (FragIndex){.rigid = false, .has_frag = true, .idx = idx_i};
    ctx->predecessors[b.aj] = b.ai;

  } else if (ctx->frag_map[b.aj].has_frag) {
    if (unlikely(ctx->frag_map[b.aj].rigid)) { return DOF_MIXED_RIGID_SEMIRIGID; }
    // Add i to j's fragment
    size_t idx_j = dofulator_get_semirigid_idx(ctx, b.aj);
    Fragment* frag = &ctx->semirigid_frags.fragments[idx_j];
    frag->n_atoms++;
    ctx->frag_map[b.ai] = (FragIndex){.rigid = false, .has_frag = true, .idx = idx_j};
    ctx->predecessors[b.ai] = b.aj;

  } else {
    // Create new fragment with i and j
    IndexResult e = fragmentlist_add_semirigid(&ctx->semirigid_frags, 2);
    if (unlikely(e.status)) { return e.status; }
    size_t frag_idx = e.idx;
    ctx->frag_map[b.ai] = (FragIndex){.rigid = false, .has_frag = true, .idx = frag_idx};
    ctx->frag_map[b.aj] = (FragIndex){.rigid = false, .has_frag = true, .idx = frag_idx};
    ctx->semirigid_frags.fragments[frag_idx].root_atom = b.ai < b.aj ? b.ai : b.aj;
    if (b.ai < b.aj) {
      ctx->predecessors[b.aj] = b.ai;
    } else {
      ctx->predecessors[b.ai] = b.aj;
    }
  }

  return DOF_SUCCESS;
}


/*******************************************************************************
 * Add atoms in bond `b` to the same rigid fragment
*/
DofulatorResult dofulator_build_rigid_fragment(Dofulator ctx, Bond b) {
  if (unlikely(b.ai < 0 || b.ai >= ctx->n_atoms || b.aj < 0 || b.aj >= ctx->n_atoms)) {
    return DOF_INDEX_ERROR;
  }
  if (unlikely(!ctx)) { return DOF_UNINITIALISED; }
  if (unlikely(b.ai == b.aj)) { return DOF_SUCCESS; }

  if (ctx->frag_map[b.ai].has_frag && ctx->frag_map[b.aj].has_frag) {
    // Both already in a fragment
    if (unlikely(!ctx->frag_map[b.ai].rigid || !ctx->frag_map[b.aj].rigid)) {
      return DOF_MIXED_RIGID_SEMIRIGID;
    }
    size_t idx_i = dofulator_get_rigid_idx(ctx, b.ai);
    size_t idx_j = dofulator_get_rigid_idx(ctx, b.aj);

    // If already in same fragment then nothing to do.
    if (idx_i == idx_j) { return DOF_SUCCESS; }

    // Two different fragments, so join them
    Fragment* frag_i = &ctx->rigid_frags.fragments[idx_i];
    Fragment* frag_j = &ctx->rigid_frags.fragments[idx_j];
    assert(frag_i->root_atom != frag_j->root_atom);
    if (frag_i->root_atom < frag_j->root_atom) {
      frag_i->n_atoms += frag_j->n_atoms;
      ctx->frag_map[b.aj].idx = idx_i;
      // Invalidate frag_j and point it to frag_i
      // Any other atoms in it will be updated later through this link
      // via `dofulator_get_rigid_idx(ctx, idx_j)`
      frag_j->invalid = true;
      frag_j->root_atom = idx_i;
    } else {
      frag_j->n_atoms += frag_i->n_atoms;
      ctx->frag_map[b.ai].idx = idx_j;
      // Invalidate frag_i and point it to frag_j
      frag_i->invalid = true;
      frag_i->root_atom = idx_j;
    }

  } else if (ctx->frag_map[b.ai].has_frag) {
    if (unlikely(!ctx->frag_map[b.ai].rigid)) { return DOF_MIXED_RIGID_SEMIRIGID; }
    // Add j to i's fragment
    size_t idx_i = dofulator_get_rigid_idx(ctx, b.ai);
    Fragment* frag = &ctx->rigid_frags.fragments[idx_i];
    frag->n_atoms++;
    if (b.aj < frag->root_atom) {
      frag->root_atom = b.aj;
    }
    ctx->frag_map[b.aj] = (FragIndex){.rigid = true, .has_frag = true, .idx = idx_i};

  } else if (ctx->frag_map[b.aj].has_frag) {
    if (unlikely(!ctx->frag_map[b.aj].rigid)) { return DOF_MIXED_RIGID_SEMIRIGID; }
    // Add i to j's fragment
    size_t idx_j = dofulator_get_rigid_idx(ctx, b.aj);
    Fragment* frag = &ctx->rigid_frags.fragments[idx_j];
    frag->n_atoms++;
    if (b.ai < frag->root_atom) {
      frag->root_atom = b.ai;
    }
    ctx->frag_map[b.ai] = (FragIndex){.rigid = true, .has_frag = true, .idx = idx_j};

  } else {
    // Create new fragment with i and j
    IndexResult e = fragmentlist_add_rigid(&ctx->rigid_frags, 2);
    if (unlikely(e.status)) { return e.status; }
    size_t frag_idx = e.idx;
    ctx->frag_map[b.ai] = (FragIndex){.rigid = true, .has_frag = true, .idx = frag_idx};
    ctx->frag_map[b.aj] = (FragIndex){.rigid = true, .has_frag = true, .idx = frag_idx};
    ctx->rigid_frags.fragments[frag_idx].root_atom = b.ai < b.aj ? b.ai : b.aj;
  }
  return DOF_SUCCESS;
}


/*******************************************************************************
 * Compare function for sorting fragments
 * Fragments with `invalid == true` are invalidated fragments which
 * have been merged into others, and hence are sorted to the end.
 * Fragments with `n_atoms == 0` are also sorted to the end to be skipped.
*/
static int fragment_sort_cmp(const void* f1, const void* f2) {
  size_t n1 = ((Fragment*)f1)->n_atoms;
  size_t n2 = ((Fragment*)f2)->n_atoms;
  if (n1 == 0 || ((Fragment*)f1)->invalid) n1 = SIZE_MAX >> 1;
  if (n2 == 0 || ((Fragment*)f2)->invalid) n2 = SIZE_MAX >> 1;
  if (n1 < n2) return -1;
  if (n1 > n2) return 1;
  return 0;
}

/*******************************************************************************
 * Consolidate merged fragments and allocate memory for buffers.
*/
DofulatorResult dofulator_finalise_fragments(Dofulator ctx) {
  // Calculate space for semi-rigid fragments
  for (size_t i = 0; i < ctx->semirigid_frags.n_fragments; ++i) {
    Fragment* frag = &ctx->semirigid_frags.fragments[i];
    if (frag->n_atoms == 0 || frag->invalid) continue;
    frag->atoms = malloc(sizeof(*frag->atoms) * frag->n_atoms);
    if (unlikely(!frag->atoms)) { return DOF_ALLOC_FAILURE; }
    fragment_set_max_modes(frag);
    if (frag->n_atoms > ctx->max_semirigid_atoms) {
      // Grow semirigid counter, zero-initialising new elements
      size_t* new_ptr = realloc(ctx->n_semirigid, sizeof(*ctx->n_semirigid) * frag->n_atoms);
      if (unlikely(!new_ptr)) { return DOF_ALLOC_FAILURE; }
      ctx->n_semirigid = new_ptr;
      memset(&ctx->n_semirigid[ctx->max_semirigid_atoms], 0,
             sizeof(*ctx->n_semirigid) * (frag->n_atoms - ctx->max_semirigid_atoms));
      ctx->max_semirigid_atoms = frag->n_atoms;
    }
    // Keep track of largest required K matrix
    if (frag->n_loops > ctx->max_K_size) {
      ctx->max_K_size = 3 * frag->n_atoms * frag->n_loops;
    }
    ctx->n_semirigid[frag->n_atoms - 1]++;
    frag->n_atoms = 0; // To be incremented back below
  }

  // Calculate space for rigid fragments
  for (size_t i = 0; i < ctx->rigid_frags.n_fragments; ++i) {
    Fragment* frag = &ctx->rigid_frags.fragments[i];
    if (frag->n_atoms == 0 || frag->invalid) continue;
    frag->atoms = malloc(sizeof(*frag->atoms) * frag->n_atoms);
    if (unlikely(!frag->atoms)) { return DOF_ALLOC_FAILURE; }
    fragment_set_max_modes(frag);
    if (frag->n_atoms > ctx->max_rigid_atoms) {
      size_t* new_ptr = realloc(ctx->n_rigid, sizeof(*ctx->n_rigid) * frag->n_atoms);
      if (unlikely(!new_ptr)) { return DOF_ALLOC_FAILURE; }
      ctx->n_rigid = new_ptr;
      memset(&ctx->n_rigid[ctx->max_rigid_atoms], 0, sizeof(*ctx->n_rigid) * (frag->n_atoms - ctx->max_rigid_atoms));
      ctx->max_rigid_atoms = frag->n_atoms;
    }
    ctx->n_rigid[frag->n_atoms - 1]++;
    frag->n_atoms = 0; // To be incremented back below
  }

  // Populate atom lists of all fragments
  // Zero out frag_idx list to track which atoms have been added to fragments
  memset(ctx->atom_frag_idx, 0, sizeof(*ctx->atom_frag_idx) * ctx->n_atoms);
  AtomTag* pred_stack = NULL;
  if (ctx->max_semirigid_atoms > 1) {
    pred_stack = malloc(sizeof(*pred_stack) * (ctx->max_semirigid_atoms - 1));
    if (unlikely(!pred_stack)) { return DOF_ALLOC_FAILURE; }
  }
  for (AtomTag i = 0; i < ctx->n_atoms; ++i) {
    FragIndex frag_idx = ctx->frag_map[i];
    if (frag_idx.has_frag) {
      Fragment* frag;
      if (frag_idx.rigid) {
        // Rigid fragment - just add the new atom
        frag = &ctx->rigid_frags.fragments[frag_idx.idx];
        ctx->atom_frag_idx[i] = frag->n_atoms;
        frag->atoms[frag->n_atoms++] = i;
      } else {

        // Semirigid fragment - need to make sure predecessors come first in the
        // fragment's atom list to simplify Jacobian building
        size_t idx = dofulator_get_semirigid_idx(ctx, i);
        frag = &ctx->semirigid_frags.fragments[idx];
        AtomTag pred = ctx->predecessors[i];
        if (i == pred) {
          // This is the root atom. Just skip if already added
          assert(i == frag->root_atom);
          if (frag->n_atoms != 0) { continue; }

        } else {
          // For non-root predecessors, atom_frag_idx will be > 0 once added
          // to a fragment. Use this to check whether predecessor is already in
          // a fragment's atom list. Add predecessors back to the root (in
          // reverse order) if needed so predecessors always come earlier in frag->atoms.
          size_t stack_idx = 0;
          while (ctx->atom_frag_idx[pred] == 0) {
            if (
              pred == frag->root_atom &&
              (frag->n_atoms > 0 ||  pred == pred_stack[stack_idx - 1])
            ) {
              // Reached the root atom and it's already added or queued, so stop
              break;
            }
            pred_stack[stack_idx++] = pred;
            pred = ctx->predecessors[pred];
          }
          // Walk back down stack and add predecessors up to i
          while (stack_idx > 0) {
            AtomTag atom = pred_stack[--stack_idx];
            ctx->atom_frag_idx[atom] = frag->n_atoms;
            frag->atoms[frag->n_atoms++] = atom;
          }
        }
        // Add atom i if it's not already in a fragment
        // Root atoms will have idx 0, but this point is only reached
        // once per root atom.
        if (ctx->atom_frag_idx[i] == 0) {
          ctx->atom_frag_idx[i] = frag->n_atoms;
          frag->atoms[frag->n_atoms++] = i;
        }
      }
    }
  }
  if (pred_stack) { free(pred_stack); }

  // Partition fragments by Jacobian size (number of atoms) for batched linalg.
  // Left-over empty fragments are placed at the end to be truncated.
  qsort(
    ctx->semirigid_frags.fragments,
    ctx->semirigid_frags.n_fragments,
    sizeof(*ctx->semirigid_frags.fragments), fragment_sort_cmp);
  qsort(
    ctx->rigid_frags.fragments,
    ctx->rigid_frags.n_fragments,
    sizeof(*ctx->rigid_frags.fragments), fragment_sort_cmp);

  // Allocate space for Jacobians

  // Semirigid
  if (ctx->max_semirigid_atoms > 0) {
    ctx->dof_buf_semirigid = malloc(sizeof(double*) * ctx->max_semirigid_atoms);
    if (unlikely(!ctx->dof_buf_semirigid)) { return DOF_ALLOC_FAILURE; }
    ctx->dof_total_buf_semirigid = malloc(sizeof(double*) * ctx->max_semirigid_atoms);
    if (unlikely(!ctx->dof_total_buf_semirigid)) { return DOF_ALLOC_FAILURE; }
    for (size_t n_atoms = 1; n_atoms <= ctx->max_semirigid_atoms; ++n_atoms) {
      size_t n_frag = ctx->n_semirigid[n_atoms-1];
      if (n_frag == 0) {
        ctx->dof_buf_semirigid[n_atoms-1] = NULL;
        ctx->dof_total_buf_semirigid[n_atoms-1] = NULL;
        continue;
      }
      ctx->dof_buf_semirigid[n_atoms-1] = malloc(sizeof(double) * n_frag * 3 * n_atoms * 3 * n_atoms);
      if (unlikely(!ctx->dof_buf_semirigid[n_atoms-1])) { return DOF_ALLOC_FAILURE; }
      ctx->dof_total_buf_semirigid[n_atoms-1] = malloc(sizeof(double) * n_frag * 3 * n_atoms);
      if (unlikely(!ctx->dof_total_buf_semirigid[n_atoms-1])) { return DOF_ALLOC_FAILURE; }
    }
  }

  // Rigid
  if (ctx->max_rigid_atoms > 0) {
    ctx->dof_buf_rigid = malloc(sizeof(double*) * ctx->max_rigid_atoms);
    if (unlikely(!ctx->dof_buf_rigid)) { return DOF_ALLOC_FAILURE; }
    ctx->dof_total_buf_rigid = malloc(sizeof(double*) * ctx->max_rigid_atoms);
    if (unlikely(!ctx->dof_total_buf_rigid)) { return DOF_ALLOC_FAILURE; }
    for (size_t n_atoms = 1; n_atoms <= ctx->max_rigid_atoms; ++n_atoms) {
      size_t n_frag = ctx->n_rigid[n_atoms-1];
      if (n_frag == 0) {
        ctx->dof_buf_rigid[n_atoms-1] = NULL;
        ctx->dof_total_buf_rigid[n_atoms-1] = NULL;
        continue;
      }
      ctx->dof_buf_rigid[n_atoms-1] = malloc(sizeof(double) * n_frag * 3 * n_atoms * 6);
      if (unlikely(!ctx->dof_buf_rigid[n_atoms-1])) { return DOF_ALLOC_FAILURE; }
      ctx->dof_total_buf_rigid[n_atoms-1] = malloc(sizeof(double) * n_frag * 3 * n_atoms);
      if (unlikely(!ctx->dof_total_buf_rigid[n_atoms-1])) { return DOF_ALLOC_FAILURE; }
    }
  }


  // Update atom to fragment mapping and store dof and dof_total pointers

  // Semirigid
  double *dof_buf, *dof_total_buf;
  size_t last_n = 0;
  for (size_t i = 0; i < ctx->semirigid_frags.n_fragments; ++i) {
    Fragment* frag = &ctx->semirigid_frags.fragments[i];
    if (frag->n_atoms == 0 || frag->invalid) {
      // Truncate left-over empty/invalid fragments from initial search
      ctx->semirigid_frags.n_fragments = i;
      break;
    }

    for (AtomTag* j = frag->atoms; j < frag->atoms + frag->n_atoms; ++j) {
      ctx->frag_map[*j].idx = i;
    }

    // Get section of dof_buf_semirigid for this fragment
    if (frag->n_atoms != last_n) {
      dof_buf = ctx->dof_buf_semirigid[frag->n_atoms-1];
      dof_total_buf = ctx->dof_total_buf_semirigid[frag->n_atoms-1];
      assert(dof_buf);
      assert(dof_total_buf);
      last_n = frag->n_atoms;
    }
    frag->dof = dof_buf;
    frag->dof_total = dof_total_buf;

    // Increment ready for next fragment
    dof_buf += (3 * frag->n_atoms * frag->n_modes);
    dof_total_buf += 3 * frag->n_atoms;
  }

  // Rigid
  last_n = 0;
  for (size_t i = 0; i < ctx->rigid_frags.n_fragments; ++i) {
    Fragment* frag = &ctx->rigid_frags.fragments[i];
    if (frag->n_atoms == 0 || frag->invalid) {
      // Truncate left-over empty fragments from initial search
      ctx->rigid_frags.n_fragments = i;
      break;
    }

    for (AtomTag* j = frag->atoms; j < frag->atoms + frag->n_atoms; ++j) {
      ctx->frag_map[*j].idx = i;
    }

    // Get section of dof_buf_rigid for this fragment
    if (frag->n_atoms != last_n) {
      dof_buf = ctx->dof_buf_rigid[frag->n_atoms-1];
      dof_total_buf = ctx->dof_total_buf_rigid[frag->n_atoms-1];
      assert(dof_buf);
      assert(dof_total_buf);
      last_n = frag->n_atoms;
    }
    frag->dof = dof_buf;
    frag->dof_total = dof_total_buf;

    // Increment ready for next fragment
    dof_buf += (3 * frag->n_atoms * frag->n_modes);
    dof_total_buf += 3 * frag->n_atoms;
  }

  // Allocate working memory for current batch size
  return dofulator_update_batch_size(ctx, ctx->batch_size);
}


/*******************************************************************************
 * Solve for DoF matrix of the given fragment.
 * TODO: batched alternative
*/
static DofulatorResult fragment_solve_dof(Dofulator ctx, Fragment* frag, double* mJ, double* Q) {
  const size_t n_atoms = frag->n_atoms;
  const size_t row_stride = frag->rigid ? 6 : 3*n_atoms;
  // Store J^T M J in Q ready for eigensolve
  cblas(dgemm)(
    CblasRowMajor, CblasTrans, CblasNoTrans,
    frag->n_modes, frag->n_modes, 3*n_atoms,
    1.0, mJ, row_stride, mJ, row_stride,
    0.0, Q, row_stride
  );

  // Find eigenvalues/eigenvectors of inertia tensor
  double work_query;
  double* Itotal = frag->dof_total; // Temporarily use dof_total for total inertia per mode
  lapack_int info = LAPACKE_dsyev_work(
    LAPACK_ROW_MAJOR, 'V', 'L',
    frag->n_modes, Q, row_stride, Itotal,
    &work_query, -1
  );
  if (unlikely(info)) { return DOF_LAPACK_ERROR; }

  if (work_query > ctx->svd_working_size) {
    ctx->svd_working_size = (lapack_int)work_query;
    double* new_ptr = realloc(ctx->svd_working, sizeof(double) * ctx->svd_working_size);
    if (unlikely(!new_ptr)) { return DOF_ALLOC_FAILURE; }
    ctx->svd_working = new_ptr;
  }
  info = LAPACKE_dsyev_work(
    LAPACK_ROW_MAJOR, 'V', 'L',
    frag->n_modes, Q, row_stride, Itotal,
    ctx->svd_working, ctx->svd_working_size
  );
  if (unlikely(info)) { return DOF_LAPACK_ERROR; }

  // Modal inertia = Q^T J^T M J Q
  // For mode m:
  //           I_m = lambda_m Q_m^T Q_m
  //               = lambda_m (= eigenvalue m)

  // Store the full transform M^1/2 J Q for later to get per-atom modal inertia
  cblas(dgemm)(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3*n_atoms, frag->n_modes, frag->n_modes,
    1.0, mJ, row_stride, Q, row_stride,
    0.0, frag->dof, row_stride
  );

  // Disregard modes with inertia less than a small fraction of the largest inertia
  // Scale by sqrt(modal inertia) so arbitrary directional DoF can be calculated from dof matrix
  double Ithresh = Itotal[frag->n_modes-1] * DBL_EPSILON;
  for (size_t m = 0; m < frag->n_modes; ++m) {
    const double Iinv = Itotal[m] > Ithresh ? 1. / sqrt(Itotal[m]) : 0.;
    cblas(dscal)(3*n_atoms, Iinv, &frag->dof[m], row_stride);
  }

  // Sum up total DoF in each axis of current reference frame.
  // Sum of 3 directions gives total DoF.
  const double* row;
  for (size_t a = 0; a < frag->n_atoms; ++a) {
    row = &frag->dof[ 3*a    * row_stride];
    frag->dof_total[3*a]   = cblas(ddot)(frag->n_modes, row, 1, row, 1);
    row += row_stride;
    frag->dof_total[3*a+1] = cblas(ddot)(frag->n_modes, row, 1, row, 1);
    row += row_stride;
    frag->dof_total[3*a+2] = cblas(ddot)(frag->n_modes, row, 1, row, 1);
  }
  return DOF_SUCCESS;
}

/*******************************************************************************
 * Find reference atoms in the fragment and store them in `frame` along with the
 * relevant unit vectors which define the reference orientation
*/
static void fragment_set_frame(const Fragment* frag, RefFrame* frame, const double x[][3], const PBC* pbc) {
  if (frag->n_atoms == 0) { return; }

  AtomTag iatom = 0;
  frame->ref_atom1 = frag->atoms[iatom];
  frame->current_rot = quat_identity();

  // Find next atom with non-zero distance from first
  frame->ref_atom2 = frag->atoms[iatom];
  double len2_a = 0., a[3];
  while (++iatom < frag->n_atoms && len2_a < 100*DBL_EPSILON) {
    vec_sub(x[frag->atoms[iatom]], x[frame->ref_atom1], a);
    pbc_wrap_shortest(pbc, a);
    len2_a = vec_dot(a, a);
  }
  --iatom;
  if (len2_a >= 100*DBL_EPSILON) {
    frame->ref_atom2 = frag->atoms[iatom];

    // Normalize vectors for simpler processing later
    const double len_a = sqrt(len2_a);
    frame->r12[0] = a[0] / len_a;
    frame->r12[1] = a[1] / len_a;
    frame->r12[2] = a[2] / len_a;
  } else {
    // Run out of atoms, so all atoms are coincident
    assert(frame->ref_atom2 == frame->ref_atom1);
    frame->ref_atom3 = frame->ref_atom2;
  }

  // 3rd atom needs non-zero cross product
  frame->ref_atom3 = frame->ref_atom2;
  double b[3], c[3];
  double len2_b = 0., len2_c = 0.;
  while (++iatom < frag->n_atoms && len2_c < 100*DBL_EPSILON) {
    vec_sub(x[frag->atoms[iatom]], x[frame->ref_atom1], b);
    pbc_wrap_shortest(pbc, b);
    len2_b = vec_dot(b, b);
    if (len2_b < 100*DBL_EPSILON) { continue; }

    vec_cross(a, b, c);
    len2_c = vec_dot(c, c);
  }
  --iatom;
  if (len2_c >= 100*DBL_EPSILON) {
    frame->ref_atom3 = frag->atoms[iatom];
    const double len_c = sqrt(len2_c);
    frame->r13_perp[0] = c[0] / len_c;
    frame->r13_perp[1] = c[1] / len_c;
    frame->r13_perp[2] = c[2] / len_c;
  } else {
    // No more available atoms, so fragment is linear
    assert(frame->ref_atom3 == frame->ref_atom2);
  }
}


/*******************************************************************************
 * Pre-calculate modal, directional DoF of rigid fragments, since these can just
 * be rotated with the fragment rather than re-calculating each time
*/
DofulatorResult dofulator_precalculate_rigid(Dofulator ctx, const double* mass, const double x[][3]) {
  if (unlikely(!ctx)) {
    return DOF_UNINITIALISED;
  }

  {
    RefFrame* new_ptr = realloc(ctx->rigid_ref_frames, sizeof(RefFrame) * ctx->rigid_frags.n_fragments);
    if (unlikely(!new_ptr)) { return DOF_ALLOC_FAILURE; }
    ctx->rigid_ref_frames = new_ptr;
  }

  for (
    size_t ifrag = 0, *n_frags = ctx->n_rigid;
    n_frags < ctx->n_rigid + ctx->max_rigid_atoms;
    ++n_frags
  ) {
    if (*n_frags == 0) { continue; }
    const size_t n_atoms = ctx->rigid_frags.fragments[ifrag].n_atoms;
    const size_t mat_size = 3*n_atoms * 6;
    double* mJ = ctx->working_buf;
    double* Q = mJ + mat_size;
    double rij[3];
    for (size_t f = 0; f < *n_frags; ++f, ++ifrag) {
      Fragment* frag = &ctx->rigid_frags.fragments[ifrag];
      double (*J)[3*n_atoms][6] = (double(*)[3*n_atoms][6])mJ;
      memset(J, 0, sizeof(*J));
      // Build fragment Jacobian
      const AtomTag root = frag->root_atom;
      for (size_t i = 0; i < frag->n_atoms; ++i) {
        const AtomTag atom = frag->atoms[i];
        const double root_mass = sqrt(mass[frag->atoms[i]]);
        (*J)[3*i  ][0] = root_mass;
        (*J)[3*i+1][1] = root_mass;
        (*J)[3*i+2][2] = root_mass;
        if (atom == root) { continue; }

        // Get vector from predecessor to atom i
        vec_sub(x[atom], x[root], rij);
        pbc_wrap_shortest(&ctx->pbc, rij);

        /* Store negative cross operator
         *  0   z   -y
         * -z   0    x
         *  y  -x    0
         */
        (*J)[3*i][4] = root_mass * rij[2];
        (*J)[3*i][5] = -root_mass * rij[1];

        (*J)[3*i+1][3] = -root_mass * rij[2];
        (*J)[3*i+1][5] = root_mass * rij[0];

        (*J)[3*i+2][3] = root_mass * rij[1];
        (*J)[3*i+2][4] = -root_mass * rij[0];
      }
      DofulatorResult e = fragment_solve_dof(ctx, frag, mJ, Q);
      if (unlikely(e)) { return e; }
      fragment_set_frame(frag, &ctx->rigid_ref_frames[ifrag], x, &ctx->pbc);
    }
  }
  return DOF_SUCCESS;
}


/*******************************************************************************
 * Get a normalized quaternion representing the rotation needed to bring the
 * fragment from its original orientation to the current one
*/
static Quaternion frame_get_rotation(const RefFrame* frame, const double x[][3], const PBC* pbc) {
  if (frame->ref_atom1 == frame->ref_atom2) {
    return quat_identity();
  }
  double a[3];
  vec_sub(x[frame->ref_atom2], x[frame->ref_atom1], a);
  pbc_wrap_shortest(pbc, a);

  // Find quaternion which rotates a onto the r12 vector
  Quaternion qa = quat_from_closest_arc(a, frame->r12);

#ifndef NDEBUG
  // Check first rotation is correct
  quat_rotate_vec(qa, a);
  vec_normalize(a);
  assert(fabs(a[0] - frame->r12[0]) < 5e-10);
  assert(fabs(a[1] - frame->r12[1]) < 5e-10);
  assert(fabs(a[2] - frame->r12[2]) < 5e-10);
#endif

  // If linear fragment, qa is all that's needed
  if (frame->ref_atom3 == frame->ref_atom2) {
    return qa;
  }

  // Apply qa rotation to b
  double b[3], c[3];
  vec_sub(x[frame->ref_atom3], x[frame->ref_atom1], b);
  pbc_wrap_shortest(pbc, b);
  quat_rotate_vec(qa, b);
  vec_cross(frame->r12, b, c);
  vec_normalize(c);

  // Get rotation around r12 axis such that c aligns with r13_perp == r12 x r13.
  // frame->r12 and r13_perp are normalized.
  // Use c x r13_perp to get rotation axis with sin(angle) scaling.
  // Use 1 + cos(angle) and then normalize to get rotation by angle/2.
  double axis_b[3];
  vec_cross(c, frame->r13_perp, axis_b);
  Quaternion qb = quat_from_vec(axis_b);
  double cos_angle = vec_dot(c, frame->r13_perp);
  qb.w = 1. + cos_angle;
  qb = quat_normalize(qb);

#ifndef NDEBUG
  // Check second rotation is correct
  quat_rotate_vec(qb, c);
  assert(fabs(c[0] - frame->r13_perp[0]) < 5e-10);
  assert(fabs(c[1] - frame->r13_perp[1]) < 5e-10);
  assert(fabs(c[2] - frame->r13_perp[2]) < 5e-10);
#endif

  // Full rotation from applying qa and then qb.
  return quat_mul(qb, qa);
}


/*******************************************************************************
 * Calculate rotation quaternion for each rigid fragment.
*/
static inline void dofulator_calculate_rigid(Dofulator ctx, const double x[][3]) {
  for (size_t ifrag = 0; ifrag < ctx->rigid_frags.n_fragments; ++ifrag) {
    ctx->rigid_ref_frames[ifrag].current_rot = frame_get_rotation(&ctx->rigid_ref_frames[ifrag], x, &ctx->pbc);
  }
}


/*******************************************************************************
 * Calculate `dof_total` and `dof` for all semirigid fragments based on given
 * masses and atom locations (indexed by atom index)
*/
static inline DofulatorResult dofulator_calculate_semirigid(Dofulator ctx, const double* mass, const double x[][3]) {
  for (
    size_t ifrag = 0, *n_frags = ctx->n_semirigid;
    n_frags < ctx->n_semirigid + ctx->max_semirigid_atoms;
    ++n_frags
  ) {
    if (*n_frags == 0) { continue; }
    const AtomTag n_atoms = ctx->semirigid_frags.fragments[ifrag].n_atoms;
    const size_t mat_size = 3*n_atoms * 3*n_atoms;
    double* mJ = ctx->working_buf;
    double* Q = mJ + MAX(mat_size, ctx->max_K_size);
    double* VT = Q + mat_size;
    double rij[3];
    for (size_t f = 0; f < *n_frags; ++f, ++ifrag) {
      Fragment* frag = &ctx->semirigid_frags.fragments[ifrag];
      double (*J)[3*n_atoms][3*n_atoms] = (double(*)[3*n_atoms][3*n_atoms]) ((frag->n_loops > 0) ? frag->dof : mJ);
      memset(J, 0, sizeof(*J));
      // Build fragment Jacobian
      for (size_t i = 0; i < frag->n_atoms; ++i) {
        const AtomTag atom = frag->atoms[i];
        const AtomTag pred = ctx->predecessors[atom];
        if (pred == atom) {
          assert(atom == frag->root_atom);
          (*J)[3*i  ][3*i  ] = 1.0;
          (*J)[3*i+1][3*i+1] = 1.0;
          (*J)[3*i+2][3*i+2] = 1.0;
        } else {
          const AtomTag pred_i = ctx->atom_frag_idx[pred];
          // Copy x,y,z rows from predecessor
          // Fragment atom list is ordered such that predecessors are always processed first.
          memcpy(&((*J)[3*i][0]), &((*J)[3*pred_i][0]), sizeof(double) * 3*n_atoms * 3);

          // Get vector from predecessor to atom i
          vec_sub(x[atom], x[pred], rij);
          pbc_wrap_shortest(&ctx->pbc, rij);

          /* Store negative cross operator
           *  0   z   -y
           * -z   0    x
           *  y  -x    0
           */
          (*J)[3*i][3*i+1] = rij[2];
          (*J)[3*i][3*i+2] = -rij[1];

          (*J)[3*i+1][3*i] = -rij[2];
          (*J)[3*i+1][3*i+2] = rij[0];

          (*J)[3*i+2][3*i] = rij[1];
          (*J)[3*i+2][3*i+1] = -rij[0];
        }
      }

      // Handle loop closures.
      // Use mJ space for K matrix, final result for null(K) stored in VT
      if (frag->n_loops > 0) {
        const lapack_int n_loops = frag->n_loops;
        const lapack_int n_modes = 3*n_atoms;
        double (*K)[frag->n_loops][n_modes] = (double(*)[frag->n_loops][n_modes])mJ;
        for (size_t k = 0; k < frag->n_loops; ++k) {
          const Bond b = frag->loop_closures[k];
          const size_t i = ctx->atom_frag_idx[b.ai];
          const size_t j = ctx->atom_frag_idx[b.aj];

          // T^T is a row vector for bonds
          vec_sub(x[b.aj], x[b.ai], rij);
          pbc_wrap_shortest(&ctx->pbc, rij);

          // Calculate J_j - J_i, store 3 x 3*n_atoms result in Q
          memcpy(Q, &((*J)[3*j][0]), sizeof(double) * n_modes*3);
          cblas(daxpy)(3*n_modes, -1.0, &((*J)[3*i][0]), 1, Q, 1);

          // Calculate row k of K = T^T (J_j - J_i)
          // Use (J_j - J_i)^T x T and store result in K[k][:]
          cblas(dgemv)(CblasRowMajor, CblasTrans, 3, n_modes, 1.0, Q, n_modes, rij, 1, 0., &((*K)[k][0]), 1);
        }

        // Find null(K) from low singular values of K
        // Temporarily store singular values in dof_total
        double* S = frag->dof_total;
        // Need to zero out S, otherwise dgesvd gives incorrect
        // results for low singular values
        memset(S, 0, sizeof(double) * n_modes);
        double work_query;
        lapack_int info = LAPACKE_dgesvd_work(
          LAPACK_ROW_MAJOR,
          'N', // Don't calculate any of U
          'A', // Calculate all of V^T
          n_loops, n_modes,
          &((*K)[0][0]), n_modes,
          S,
          NULL, frag->n_loops,
          VT, n_modes,
          &work_query,
          -1
        );
        if (unlikely(info)) { return DOF_LAPACK_ERROR; }

        // Need more working space, so allocate it
        if (work_query > ctx->svd_working_size) {
          ctx->svd_working_size = (lapack_int)work_query;
          double* new_ptr = realloc(ctx->svd_working, sizeof(double) * ctx->svd_working_size);
          if (unlikely(!new_ptr)) { return DOF_ALLOC_FAILURE; }
          ctx->svd_working = new_ptr;
        }
        info = LAPACKE_dgesvd_work(
          LAPACK_ROW_MAJOR,
          'N', // Don't calculate any of U
          'A', // Calculate all of V^T
          n_loops, n_modes,
          &((*K)[0][0]), n_modes,
          S,
          NULL, frag->n_loops,
          VT, n_modes,
          ctx->svd_working,
          ctx->svd_working_size
        );
        if (unlikely(info)) { return DOF_LAPACK_ERROR; }

        // Find rank of null(K)
        // Singular values in descending order
        // Take max with a small cut-off in case all loop closures
        // are redundant constraints
        double thresh;
        if (likely(ctx->null_space_thresh == 0.0)) {
          thresh = DBL_EPSILON * MAX(n_modes, n_loops);
        } else {
          thresh = ctx->null_space_thresh;
        }
        thresh = MAX(thresh * S[0], DBL_EPSILON);
        size_t rank_K = 0;
        while (rank_K < (size_t)n_modes && S[rank_K] > thresh) { ++rank_K; }

        // Project Jacobian onto null(K).
        // null(K) = ( VT[rank_K:3*n_atoms][..] )^T
        cblas(dgemm)(
          CblasRowMajor, CblasNoTrans, CblasTrans,
          n_modes, n_modes - rank_K, n_modes,
          1.0, &((*J)[0][0]), n_modes,
          &VT[rank_K * n_modes], n_modes,
          0.0, mJ, n_modes
        );

        frag->n_modes = n_modes - rank_K;
      }

      // Scale rows by sqrt(mass) of corresponding atom
      for (AtomTag a = 0; a < n_atoms; ++a) {
        // TODO: pre-calculate sqrt(mass) and store on fragment?
        double root_mass = sqrt(mass[frag->atoms[a]]);
        cblas(dscal)(frag->n_modes, root_mass, &mJ[3*a*3*n_atoms], 1);
        cblas(dscal)(frag->n_modes, root_mass, &mJ[(3*a+1)*3*n_atoms], 1);
        cblas(dscal)(frag->n_modes, root_mass, &mJ[(3*a+2)*3*n_atoms], 1);
      }

      DofulatorResult e = fragment_solve_dof(ctx, frag, mJ, Q);
      if (unlikely(e)) { return e; }
    }
  }

  return DOF_SUCCESS;

  // PLANNING:
  // For each semirigid fragment in batch, need:
  //  - Jacobian space (mJ): (3*max_semirigid_atoms)^2
  //  - Eigenvector space:
  //    + Q: (3*max_semirigid_atoms)^2
  //    + dof_total: 3*max_semirigid_atoms
  //  - Auxilliary space for n_closures loop closures:
  //    + space for K matrix construction:
  //      = T: 3 <- Use same space as VT since can be overwritten
  //      = (J_j - J_i): 3*(3*max_semirigid_atoms)
  //        * Use same space as VT since can be overwritten
  //        * Probably fastest to use daxpy on full batch to do subtraction
  //          and store in VT memory, then dgemv to fill K?
  //      = K: n_closures * 3*max_semirigid_atoms
  //        * Fill by re-using T, (J_j - Ji) memory space for each closure
  //    + space for SVD of K:
  //      = U: none - pass NULL for U arg
  //      = S: n_closures
  //      = VT: (3*max_semirigid_atoms)^2
  //      = WORK: get from dgesvd call
  //    + space for J * V[:, rank(S):]
  //      = Can store in eigenvector space and pointer swap back into J
  //  (use eigenvector space since this comes first?)
}


/*******************************************************************************
 * Calculate DoF of all fragments in the context
*/
DofulatorResult dofulator_calculate(Dofulator ctx, const double* mass, const double x[][3]) {
  dofulator_calculate_rigid(ctx, x);
  return dofulator_calculate_semirigid(ctx, mass, x);
}


/*******************************************************************************
 * Get the total DoF of the atom with index `atom_idx`
*/
double dofulator_get_dof_atom(const struct Dofulator* ctx, AtomTag atom_idx) {
  if (atom_idx >= ctx->n_atoms) {
    return 0.;
  }
  FragIndex frag_idx = ctx->frag_map[atom_idx];
  if (!frag_idx.has_frag) {
    // TODO: Account for global CoM velocity constraint
    return 3.;
  }

  const size_t i = ctx->atom_frag_idx[atom_idx];
  Fragment* frag;
  if (frag_idx.rigid) {
    frag = &ctx->rigid_frags.fragments[frag_idx.idx];
  } else {
    frag = &ctx->semirigid_frags.fragments[frag_idx.idx];
  }
  return frag->dof_total[3*i] + frag->dof_total[3*i + 1] + frag->dof_total[3*i + 2];
}


/*******************************************************************************
 * Get the directional DoF of the atom with index `atom_idx`
*/
void dofulator_get_dof_atom_directional(const struct Dofulator* restrict ctx, AtomTag atom_idx, double dof[restrict 3]) {
  dof[0] = dof[1] = dof[2] = 0.;
  if (atom_idx >= ctx->n_atoms) {
    return;
  }
  FragIndex frag_idx = ctx->frag_map[atom_idx];
  if (!frag_idx.has_frag) {
    // TODO: Account for global CoM velocity constraint
    dof[0] = dof[1] = dof[2] = 1.;
    return;
  }
  Fragment* frag;
  if (frag_idx.rigid) {
    frag = &ctx->rigid_frags.fragments[frag_idx.idx];
  } else {
    frag = &ctx->semirigid_frags.fragments[frag_idx.idx];
  }
  const size_t i = ctx->atom_frag_idx[atom_idx];
  if (frag_idx.rigid) {
    // Rotate into current frame for rigid fragments
    double x[3] = {1., 0., 0.};
    double y[3] = {0., 1., 0.};
    double z[3] = {0., 0., 1.};
    Quaternion q = ctx->rigid_ref_frames[frag_idx.idx].current_rot;
    quat_rotate_vec(q, x);
    quat_rotate_vec(q, y);
    quat_rotate_vec(q, z);

    // Sum up total DoF in each direction
    // For direction e_a, DoF = sum_{m in modes} e_m (mJQ_i)^T e_a^T e_a mJQ_i e_m / lambda_m
    // where mJQ_i = rows of sqrt(M)JQ corresponding to atom i, and lambda_m = modal inertia.
    // frag->dof matrix stores sqrt(m)JQ/sqrt(lambda_m)
    double d, rdx, rdy, rdz;
    for (size_t m = 0; m < 6; ++m) {
      rdx = frag->dof[ 3*i    * 6 + m];
      rdy = frag->dof[(3*i+1) * 6 + m];
      rdz = frag->dof[(3*i+2) * 6 + m];
      d = x[0]*rdx + x[1]*rdy + x[2]*rdz;
      dof[0] += d*d;
      d = y[0]*rdx + y[1]*rdy + y[2]*rdz;
      dof[1] += d*d;
      d = z[0]*rdx + z[1]*rdy + z[2]*rdz;
      dof[2] += d*d;
    }
  } else {
    dof[0] = frag->dof_total[3*i];
    dof[1] = frag->dof_total[3*i+1];
    dof[2] = frag->dof_total[3*i+2];
  }

}


/*******************************************************************************
 * Get the total DoF of the `n_atoms` atoms with indices in `atoms`.
 * Allocates `dof` if NULL is passed in.
 * Returns NULL if `dof` is NULL and allocation fails.
*/
double* dofulator_get_dof_atoms(
  const struct Dofulator* ctx,
  const size_t n_atoms,
  const AtomTag* atoms,
  double* restrict dof
) {
  if (dof == NULL) {
    dof = malloc(sizeof(double) * n_atoms);
    if (!dof) return NULL;
  }
  if (atoms) {
    for (size_t i = 0; i < n_atoms; ++i) {
      dof[i] = dofulator_get_dof_atom(ctx, atoms[i]);
    }
  } else {
    for (size_t i = 0; i < n_atoms; ++i) {
      dof[i] = dofulator_get_dof_atom(ctx, i);
    }
  }
  return dof;
}


/*******************************************************************************
 * Get the directional DoF of the `n_atoms` atoms with indices in `atoms`.
 * Allocates `dof` if NULL is passed in.
 * Returns NULL if `dof` is NULL and allocation fails.
*/
double* dofulator_get_dof_atoms_directional(
  const struct Dofulator* restrict ctx,
  const size_t n_atoms,
  const AtomTag* restrict atoms,
  double dof[restrict][3]
) {
  if (dof == NULL) {
    dof = (double(*)[3])malloc(sizeof(double) * n_atoms * 3);
    if (!dof) return NULL;
  }
  if (atoms) {
    for (size_t i = 0; i < n_atoms; ++i) {
      dofulator_get_dof_atom_directional(ctx, atoms[i], dof[i]);
    }
  } else {
    for (size_t i = 0; i < n_atoms; ++i) {
      dofulator_get_dof_atom_directional(ctx, i, dof[i]);
    }
  }
  return &dof[0][0];
}


/*******************************************************************************
 * Get an iterator over the list of rigid bodies
*/
FragmentListIter dofulator_get_rigid_fragments(const struct Dofulator* ctx) {
  return (FragmentListIter){
    .fragments = ctx->rigid_frags.fragments,
    .n_fragments = ctx->rigid_frags.n_fragments,
    .idx = 0,
  };
}

/*******************************************************************************
 * Get an iterator over the list of semi-rigid fragments
*/
FragmentListIter dofulator_get_semirigid_fragments(const struct Dofulator* ctx) {
  return (FragmentListIter){
    .fragments = ctx->semirigid_frags.fragments,
    .n_fragments = ctx->semirigid_frags.n_fragments,
    .idx = 0,
  };
}


/*******************************************************************************
 * Get the next fragment from the iterator
*/
const Fragment* fragmentlist_iter_next(FragmentListIter* iter) {
  if (iter->idx < iter->n_fragments) {
    return &iter->fragments[iter->idx++];
  } else {
    return NULL;
  }
}


/*******************************************************************************
 * Get a read-only view of atoms in a given fragment
*/
AtomListView fragment_get_atoms(const Fragment* frag) {
  return (AtomListView){.atoms = frag->atoms, .n_atoms = frag->n_atoms};
}

