#include <assert.h>
#include <float.h>
#include <math.h>
#include <openblas/cblas.h>
#include <openblas/lapacke.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "cblas_lapacke.h"
#include "dofulator.h"
#include "quaternion.h"
#include "vec3.h"

typedef struct Fragment {
  uint32_t n_atoms;         // Number of atoms in the fragment
  uint32_t n_modes;         // Maximum number of modes this fragment could have
  uint32_t rigid : 1;       // 1 = rigid, 0 = semirigid
  uint32_t invalid : 1;     // 1 = fragment has been invalidated, in which case root_atom
                            // is the index of the fragment it was merged into
  uint32_t n_loops : 30;    // Number of kinematic loops (length of `loop_closures`)

  size_t root_atom;         // Root atom of the fragment's tree.
                            // For merged fragments, this is the idx of the parent
                            // fragment in the context's `frag_map`.

  AtomTag* atoms;           // Indices of atoms in this fragment. Stored in ascending order.
  Bond* loop_closures;      // Rigid bonds that create kinematic loops

  double* dof;              // 3*n_atoms x n_modes matrix of M^1/2 J Q
  double* dof_total;        // 3*n_atoms vector of total direciontal dof
} Fragment;


// Reference frame for rigid fragments
typedef struct RefFrame {
  size_t ref_atom1, ref_atom2, ref_atom3; // Atoms to use to construct the reference frame for rigid fragments
                                          // ref_atom2 == ref_atom3 if only 2 atoms
  double r12[3], r13_perp[3]; // r12 and (r12 x r13) vectors in the original coordinates
  Quaternion current_rot;
} RefFrame;

typedef struct FragIndex {
  size_t has_frag : 1;
  size_t rigid : 1;
  size_t idx : (sizeof(size_t) - 2);
} FragIndex;

typedef struct FragmentList {
  size_t capacity;
  size_t n_fragments;   // Number of fragments
  Fragment* fragments;  // Individual fragments
} FragmentList;


struct Dofulator {
  AtomTag n_atoms;                // Number of atoms in frag_map
  FragIndex* frag_map;            // Index of fragment corresponding to given atom
  AtomTag* predecessors;          // Predecessor for each atom (self for fragment root nodes)
  size_t* atom_frag_idx;          // Index of the atom within the fragment. Gives corresponding rows in Jacobian
  size_t* n_semirigid;            // Number of semirigid fragments indexed by (number of atoms - 1)
  size_t* n_rigid;                // Number of rigid fragments indexed by number of atoms
  size_t max_semirigid_atoms;     // Maximum number of atoms in any semirigid fragment
  size_t max_rigid_atoms;         // Maximum number of atoms in any rigid fragment
  size_t batch_size;              // Minimum number of fragments to process at once.
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
  double null_space_thresh;       // Fraction of the maximum singular value below which singular values are treated as 0
  lapack_int svd_working_size;    // Size of SVD/eigenvector working memory
};

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

static inline Fragment fragment_create(AtomTag n_atoms, bool rigid) {
  return (Fragment){
    .n_atoms = n_atoms,
    .rigid = rigid,
    .n_loops = 0,
  };
}

void fragment_set_max_modes(Fragment* frag) {
  frag->n_modes = frag->rigid ? 6 : 3 * frag->n_atoms;
}

// Add a new rigid fragment with `n_atoms` atoms to `list`.
#define fragmentlist_add_rigid(list, n_atoms) fragmentlist_add_new((list), (n_atoms), true)

// Add a new semirigid fragment with `n_atoms` atoms to `list`.
#define fragmentlist_add_semirigid(list, n_atoms) fragmentlist_add_new((list), (n_atoms), false)

// Add a new fragment with `n_atoms` atoms to `list`.
// `rigid` should be `true` for rigid fragments, or `false` for semirigid.
size_t fragmentlist_add_new(FragmentList* list, size_t n_atoms, bool rigid) {
  // Find the next empty slot
  // Allocate if needed
  if (list->n_fragments >= list->capacity) {
    // List needs to grow, so pick a new capacity and extend it
    list->capacity = list->capacity < 4 ? 4 : list->capacity * 2;
    list->fragments = realloc(list->fragments, sizeof(Fragment) * list->capacity);
  }

  // Create the new fragment.
  // Allocation and filling of atom list/buffers comes later.
  size_t idx = list->n_fragments++;
  list->fragments[idx] = fragment_create(n_atoms, rigid);

  return idx;
}

void fragmentlist_destroy(FragmentList* self) {
  if (!self) return;
  if (self->fragments) free(self->fragments);
  self->n_fragments = 0;
  self->capacity = 0;
}

// Get index of the top-level fragment containing fragment `idx`.
// Used to get updated fragment index from an invalidated (merged) fragment.
size_t fragmentlist_get_idx(FragmentList fraglist, size_t idx) {
  if (!fraglist.fragments[idx].invalid) {
    return idx;
  }
  return fragmentlist_get_idx(fraglist, fraglist.fragments[idx].root_atom);
}

// Create a dofulator context with capacity for `n_atoms` total atoms (could
// contain many rigid or semi-rigid fragments)
Dofulator dofulator_create(AtomTag n_atoms) {
  struct Dofulator ctx = {
    .n_atoms = n_atoms,
    .batch_size = DOFULATOR_DEFAULT_BATCH_SIZE,
    .frag_map = (FragIndex*)calloc(n_atoms, sizeof(FragIndex)),
    .predecessors = (AtomTag*)malloc(sizeof(AtomTag) * n_atoms),
    .atom_frag_idx = (size_t*)malloc(sizeof(size_t) * n_atoms),
    .null_space_thresh = 0.001,   // This seems quite high, but helps with stiff systems
  };
  for (AtomTag i = 0; i < n_atoms; ++i) {
    ctx.predecessors[i] = i;
  }
  Dofulator out = malloc(sizeof(ctx));
  // TODO: check for/handle errors
  *out = ctx;
  return out;
}

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

// Set the batch size and reallocate working memory
// Working memory:
// Max. Jacobian size = MAX((max_semirigid_atoms * 3)^2, (max_rigid_atoms*3 * 6))
// Max. Eigenvector size = MAX((max_semirigid_atoms * 3)^2, 6^2)
//
// Solve process:
// Build J (in mJQ for loops, or in mJ for trees/rigid)
// (loops only)
//    build K from J (can use mJ scratch space for result, VT scratch for J_j - J_i, Q scratch for T^T) -> batched dgemm for each T^T (J_j - J_i)
//    find null(K) (store in VT scratch space, singular values in dof_total)
//      -> dgesvd in loop over batch. Pass NULL for U.
//      -> Also needs WORK storage (get size from dgesvd call, can allocate for single fragment and re-use)
//    project J onto null(K) (store in mJ scratch space) -> batched dgemm
// Calculate sqrt(m) J (directly in scratch mJ space) -> dscal (batched?) modifies in-place
// Calculate mJ^T mJ (store in Q scratch space) -> batched dgemm
// Calculate eigenvectors (result replaces Q) -> dsyev in loop over batch, eigenvalues in Itotal
// Calculate mJQ = mJ (mJ scratch space) * Q (Q scratch space) -> batched dgemm
// Calculate dof[i][m] = mJQ[i][m] mJQ[i][m] / Itotal[m]
//
// Storage requirements per batch:
// - mJ space (Jacobian size)
// - Q  space (Eigenvector size)
// - VT space (Jacobian size - semi-rigid only)
// - WORK space for SVD (pick largest required by any loop-closing fragments)
void dofulator_update_batch_size(Dofulator ctx, size_t batch_size) {
  ctx->batch_size = batch_size;
  size_t semirigid_size = 3 * ctx->max_semirigid_atoms * 3 * ctx->max_semirigid_atoms;
  size_t max_jacobian = MAX( 3 * ctx->max_rigid_atoms * 6, semirigid_size);
  size_t max_eigenvector = MAX(36, semirigid_size);
  size_t max_scratch_size_single = max_jacobian + max_eigenvector + semirigid_size;
  ctx->working_buf_size = batch_size * max_scratch_size_single;
  ctx->working_buf = realloc(ctx->working_buf, sizeof(double) * ctx->working_buf_size);
  // TODO: better error handling
  assert(ctx->working_buf);
}

// Get the index of the semirigid fragment which contains the atom with indiex `atom_idx`.
// May mutate `ctx->frag_map` if it points to an invalidated (merged) fragment.
size_t dofulator_get_semirigid_idx(Dofulator ctx, size_t atom_idx) {
  size_t guess = ctx->frag_map[atom_idx].idx;
  size_t idx = fragmentlist_get_idx(ctx->semirigid_frags, guess);
  if (idx != guess) {
    ctx->frag_map[atom_idx].idx = idx;
  }
  return idx;
}

// Get the index of the rigid fragment which contains the atom with indiex `atom_idx`.
// May mutate `ctx->frag_map` if it points to an invalidated (merged) fragment.
size_t dofulator_get_rigid_idx(Dofulator ctx, size_t atom_idx) {
  size_t guess = ctx->frag_map[atom_idx].idx;
  size_t idx = fragmentlist_get_idx(ctx->rigid_frags, guess);
  if (idx != guess) {
    ctx->frag_map[atom_idx].idx = idx;
  }
  return idx;
}

// Update list of semirigid fragments to account for a rigid bond.
void dofulator_add_rigid_bond(Dofulator ctx, Bond b) {
  assert(b.ai < ctx->n_atoms && b.aj < ctx->n_atoms);
  assert(ctx);
  if (b.ai == b.aj) { return; }

  if (ctx->frag_map[b.ai].has_frag && ctx->frag_map[b.aj].has_frag) {
    assert(!ctx->frag_map[b.ai].rigid);
    assert(!ctx->frag_map[b.aj].rigid);
    size_t idx_i = dofulator_get_semirigid_idx(ctx, b.ai);
    size_t idx_j = dofulator_get_semirigid_idx(ctx, b.aj);
    // Both already in a fragment
    if (idx_i == idx_j) {
      // Already in same fragment, so this is a loop closure
      Fragment* frag = &ctx->semirigid_frags.fragments[idx_i];
      frag->n_loops++;
      frag->loop_closures = realloc(frag->loop_closures, sizeof(Bond) * frag->n_loops);
      // TODO: better error handling
      assert(frag->loop_closures);
      frag->loop_closures[frag->n_loops - 1] = b;
      return; // Return here so predecessors aren't changed
    } else {
      // Two different fragments, so join them
      Fragment* frag_i = &ctx->semirigid_frags.fragments[idx_i];
      Fragment* frag_j = &ctx->semirigid_frags.fragments[idx_j];
      assert(frag_i->root_atom != frag_j->root_atom);
      if (frag_i->root_atom < frag_j->root_atom) {
        frag_i->n_atoms += frag_j->n_atoms;
        ctx->frag_map[b.aj].idx = idx_i;
        // Invalidate frag_j and point it to frag_i
        // Any other atoms in it will be updated later through this link
        // via `dofulator_get_semirigid_idx(ctx, idx_j)`
        frag_j->invalid = true;
        frag_j->root_atom = idx_i;
      } else {
        frag_j->n_atoms += frag_i->n_atoms;
        ctx->frag_map[b.ai].idx = idx_j;
        // Invalidate frag_i and point it to frag_j
        frag_i->invalid = true;
        frag_i->root_atom = idx_j;
      }
    }

  } else if (ctx->frag_map[b.ai].has_frag) {
    assert(!ctx->frag_map[b.ai].rigid);
    // Add j to i's fragment
    size_t idx_i = dofulator_get_semirigid_idx(ctx, b.ai);
    Fragment* frag = &ctx->semirigid_frags.fragments[idx_i];
    frag->n_atoms++;
    if (b.aj < frag->root_atom) {
      frag->root_atom = b.aj;
    }
    ctx->frag_map[b.aj] = (FragIndex){.rigid = false, .has_frag = true, .idx = idx_i};
    
  } else if (ctx->frag_map[b.aj].has_frag) {
    assert(!ctx->frag_map[b.aj].rigid);
    // Add i to j's fragment
    size_t idx_j = dofulator_get_semirigid_idx(ctx, b.aj);
    Fragment* frag = &ctx->semirigid_frags.fragments[idx_j];
    frag->n_atoms++;
    if (b.ai < frag->root_atom) {
      frag->root_atom = b.ai;
    }
    ctx->frag_map[b.ai] = (FragIndex){.rigid = false, .has_frag = true, .idx = idx_j};

  } else {
    // Create new fragment with i and j
    size_t frag_idx = fragmentlist_add_semirigid(&ctx->semirigid_frags, 2);
    ctx->frag_map[b.ai] = (FragIndex){.rigid = false, .has_frag = true, .idx = frag_idx};
    ctx->frag_map[b.aj] = (FragIndex){.rigid = false, .has_frag = true, .idx = frag_idx};
    ctx->semirigid_frags.fragments[frag_idx].root_atom = b.ai < b.aj ? b.ai : b.aj;
  }

  if (b.ai < b.aj) {
    ctx->predecessors[b.aj] = b.ai;
  } else {
    ctx->predecessors[b.ai] = b.aj;
  }
}

// Add atoms in bond `b` to a rigid fragment
void dofulator_build_rigid_fragment(Dofulator ctx, Bond b) {
  assert(b.ai < ctx->n_atoms && b.aj < ctx->n_atoms);
  assert(ctx);
  if (b.ai == b.aj) { return; }

  if (ctx->frag_map[b.ai].has_frag && ctx->frag_map[b.aj].has_frag) {
    assert(ctx->frag_map[b.ai].rigid);
    assert(ctx->frag_map[b.aj].rigid);
    size_t idx_i = dofulator_get_rigid_idx(ctx, b.ai);
    size_t idx_j = dofulator_get_rigid_idx(ctx, b.aj);
    // Both already in a fragment
    if (idx_i == idx_j) {
      // Already in same fragment, so nothing to do.
      return;
    } else {
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
    }

  } else if (ctx->frag_map[b.ai].has_frag) {
    assert(ctx->frag_map[b.ai].rigid);
    // Add j to i's fragment
    size_t idx_i = dofulator_get_rigid_idx(ctx, b.ai);
    Fragment* frag = &ctx->rigid_frags.fragments[idx_i];
    frag->n_atoms++;
    if (b.aj < frag->root_atom) {
      frag->root_atom = b.aj;
    }
    ctx->frag_map[b.aj] = (FragIndex){.rigid = true, .has_frag = true, .idx = idx_i};
    
  } else if (ctx->frag_map[b.aj].has_frag) {
    assert(ctx->frag_map[b.aj].rigid);
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
    size_t frag_idx = fragmentlist_add_rigid(&ctx->rigid_frags, 2);
    ctx->frag_map[b.ai] = (FragIndex){.rigid = true, .has_frag = true, .idx = frag_idx};
    ctx->frag_map[b.aj] = (FragIndex){.rigid = true, .has_frag = true, .idx = frag_idx};
    ctx->rigid_frags.fragments[frag_idx].root_atom = b.ai < b.aj ? b.ai : b.aj;
  }
}

// Compare function for sorting fragments
// Fragments with `invalid == true` are invalidated fragments which
// have been merged into others, and hence are sorted to the end.
// Fragments with `n_atoms == 0` are also sorted to the end to be skipped.
int fragment_sort_cmp(const void* f1, const void* f2) {
  size_t n1 = ((Fragment*)f1)->n_atoms;
  size_t n2 = ((Fragment*)f2)->n_atoms;
  if (n1 == 0 || ((Fragment*)f1)->invalid) n1 = SIZE_MAX >> 1;
  if (n2 == 0 || ((Fragment*)f2)->invalid) n2 = SIZE_MAX >> 1;
  if (n1 < n2) return -1;
  if (n1 > n2) return 1;
  return 0;
}

// Finalise fragments by consolidating merged fragments
// and allocating memory for buffers.
void dofulator_finalise_fragments(Dofulator ctx) {
  // Calculate space for semi-rigid fragments
  for (size_t i = 0; i < ctx->semirigid_frags.n_fragments; ++i) {
    Fragment* frag = &ctx->semirigid_frags.fragments[i];
    if (frag->n_atoms == 0 || frag->invalid) continue;
    frag->atoms = malloc(sizeof(*frag->atoms) * frag->n_atoms);
    fragment_set_max_modes(frag);
    if (frag->n_atoms > ctx->max_semirigid_atoms) {
      // Grow semirigid counter, zero-initialising new elements
      ctx->n_semirigid = realloc(ctx->n_semirigid, sizeof(*ctx->n_semirigid) * frag->n_atoms);
      // TODO: better error handling
      assert(ctx->n_semirigid);
      memset(&ctx->n_semirigid[ctx->max_semirigid_atoms], 0,
             sizeof(*ctx->n_semirigid) * (frag->n_atoms - ctx->max_semirigid_atoms));
      ctx->max_semirigid_atoms = frag->n_atoms;
    }
    ctx->n_semirigid[frag->n_atoms - 1]++;
    frag->n_atoms = 0; // To be incremented back below
  }

  // Calculate space for rigid fragments
  for (size_t i = 0; i < ctx->rigid_frags.n_fragments; ++i) {
    Fragment* frag = &ctx->rigid_frags.fragments[i];
    if (frag->n_atoms == 0 || frag->invalid) continue;
    frag->atoms = malloc(sizeof(*frag->atoms) * frag->n_atoms);
    fragment_set_max_modes(frag);
    if (frag->n_atoms > ctx->max_rigid_atoms) {
      ctx->n_rigid = realloc(ctx->n_rigid, sizeof(*ctx->n_rigid) * frag->n_atoms);
      // TODO: better error handling
      assert(ctx->n_rigid);
      memset(&ctx->n_rigid[ctx->max_rigid_atoms], 0, sizeof(*ctx->n_rigid) * (frag->n_atoms - ctx->max_rigid_atoms));
      ctx->max_rigid_atoms = frag->n_atoms;
    }
    ctx->n_rigid[frag->n_atoms - 1]++;
    frag->n_atoms = 0; // To be incremented back below
  }

  // Populate atom lists of all fragments
  for (AtomTag i = 0; i < ctx->n_atoms; ++i) {
    FragIndex frag_idx = ctx->frag_map[i];
    if (frag_idx.has_frag) {
      Fragment* frag;
      if (frag_idx.rigid) {
        frag = &ctx->rigid_frags.fragments[frag_idx.idx];
      } else {
        size_t idx = dofulator_get_semirigid_idx(ctx, i);
        frag = &ctx->semirigid_frags.fragments[idx];
      }
      ctx->atom_frag_idx[i] = frag->n_atoms;
      frag->atoms[frag->n_atoms++] = i;
    }
  }

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
  ctx->dof_buf_semirigid = malloc(sizeof(double*) * ctx->max_semirigid_atoms);
  ctx->dof_total_buf_semirigid = malloc(sizeof(double*) * ctx->max_semirigid_atoms);
  for (size_t n_atoms = 1; n_atoms <= ctx->max_semirigid_atoms; ++n_atoms) {
    size_t n_frag = ctx->n_semirigid[n_atoms-1];
    if (n_frag == 0) {
      ctx->dof_buf_semirigid[n_atoms-1] = NULL;
      ctx->dof_total_buf_semirigid[n_atoms-1] = NULL;
      continue;
    }
    ctx->dof_buf_semirigid[n_atoms-1] = malloc(sizeof(double) * n_frag * 3 * n_atoms * 3 * n_atoms);
    ctx->dof_total_buf_semirigid[n_atoms-1] = malloc(sizeof(double) * n_frag * 3 * n_atoms);
  }
  
  // Rigid
  ctx->dof_buf_rigid = malloc(sizeof(double*) * ctx->max_rigid_atoms);
  ctx->dof_total_buf_rigid = malloc(sizeof(double*) * ctx->max_rigid_atoms);
  for (size_t n_atoms = 1; n_atoms <= ctx->max_rigid_atoms; ++n_atoms) {
    size_t n_frag = ctx->n_rigid[n_atoms-1];
    if (n_frag == 0) {
      ctx->dof_buf_rigid[n_atoms-1] = NULL;
      ctx->dof_total_buf_rigid[n_atoms-1] = NULL;
      continue;
    }
    ctx->dof_buf_rigid[n_atoms-1] = malloc(sizeof(double) * n_frag * 3 * n_atoms * 6);
    ctx->dof_total_buf_rigid[n_atoms-1] = malloc(sizeof(double) * n_frag * 3 * n_atoms);
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
  dofulator_update_batch_size(ctx, ctx->batch_size);
}


// Solve for DoF matrix of the given fragment.
// TODO: batched alternative
void fragment_solve_dof(Dofulator ctx, Fragment* frag, double* mJ, double* Q) {
  const size_t n_atoms = frag->n_atoms;
  const size_t row_stride = frag->rigid ? 6 : 3*n_atoms;
  // Store J^T M J in Q ready for eigensolve
  cblas_dgemm(
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
  // TODO: better error handling
  assert(info == 0);
  if (work_query > ctx->svd_working_size) {
    ctx->svd_working_size = (lapack_int)work_query;
    ctx->svd_working = realloc(ctx->svd_working, sizeof(double) * ctx->svd_working_size);
    // TODO: better error handling
    assert(ctx->svd_working);
  }
  info = LAPACKE_dsyev_work(
    LAPACK_ROW_MAJOR, 'V', 'L',
    frag->n_modes, Q, row_stride, Itotal,
    ctx->svd_working, ctx->svd_working_size
  );
  // TODO: better error handling
  assert(info == 0);

  // Modal inertia = Q^T J^T M J Q
  // For mode m:
  //           I_m = lambda_m Q_m^T Q_m
  //               = lambda_m (eigenvalue m)

  // Store the full transform M^1/2 J Q for later to get per-atom modal inertia
  cblas_dgemm(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3*n_atoms, frag->n_modes, frag->n_modes,
    1.0, mJ, row_stride, Q, row_stride,
    0.0, frag->dof, row_stride
  );

  // Calculate directional DoF per atom per mode
  for (size_t r = 0; r < 3*n_atoms; ++r) {
    for (size_t m = 0; m < frag->n_modes; ++m) {
      // TODO: blas routine for this?
      frag->dof[r*row_stride + m] *= frag->dof[r*row_stride + m];
    }
  }
  for (size_t m = 0; m < frag->n_modes; ++m) {
    const double Iinv = Itotal[m] > 100. * DBL_EPSILON ? 1. / Itotal[m] : 0.;
    cblas_dscal(3*n_atoms, Iinv, &frag->dof[m], row_stride);
  }

  // Sum up total DoF in each direction
  for (size_t a = 0; a < frag->n_atoms; ++a) {
    frag->dof_total[3*a]   = cblas_dsum(frag->n_modes, &frag->dof[ 3*a    * row_stride], 1);
    frag->dof_total[3*a+1] = cblas_dsum(frag->n_modes, &frag->dof[(3*a+1) * row_stride], 1);
    frag->dof_total[3*a+2] = cblas_dsum(frag->n_modes, &frag->dof[(3*a+2) * row_stride], 1);
  }
}

// Get a normalized quaternion representing the rotation needed to bring the fragment
// from its original orientation to the current one
Quaternion frame_get_rotation(const RefFrame* frame, const double (*x)[3]) {
  if (frame->ref_atom1 == frame->ref_atom2) {
    return quat_identity();
  }
  double a[3];
  vec_sub(x[frame->ref_atom2], x[frame->ref_atom1], a);

  // Find quaternion which rotates a onto the r12 vector
  Quaternion qa = quat_from_closest_arc(a, frame->r12);

  // If linear fragment, qa is all that's needed
  if (frame->ref_atom3 == frame->ref_atom2) {
    return qa;
  }

  // Apply qa rotation to b
  double b[3], c[3];
  vec_sub(x[frame->ref_atom3], x[frame->ref_atom1], b);
  quat_rotate_vec(qa, b);
  vec_cross(frame->r12, b, c);
  vec_normalize(c);

  // Get rotation around r12 axis such that c aligns with r13_perp == r12 x r13.
  // frame->r12 and r13_perp are normalized.
  Quaternion qb = quat_from_vec(frame->r12);
  double cos_angle = vec_dot(c, frame->r13_perp);
  if (fabs(1. + cos_angle) > 100.*DBL_EPSILON) {
    // Need non-zero axis if w is zero or normalizing gives NaN
    double sin_angle = sqrt(1 - cos_angle*cos_angle);
    qb.x *= sin_angle;
    qb.y *= sin_angle;
    qb.z *= sin_angle;
  }
  qb.w = 1. + cos_angle;
  qb = quat_normalize(qb);

  return quat_mul(qb, qa);
}

// Find reference atoms in the fragment and store them in `frame` along with the
// relevant unit vectors which define the reference orientation
void fragment_set_frame(const Fragment* frag, RefFrame* frame, const double (*x)[3]) {
  if (frag->n_atoms == 0) { return; }

  AtomTag iatom = 0;
  frame->ref_atom1 = frag->atoms[iatom];
  frame->current_rot = quat_identity();

  // Find next atom with non-zero distance from first
  frame->ref_atom2 = frag->atoms[iatom];
  double len2_a = 0., a[3];
  while (++iatom < frag->n_atoms && len2_a < 100*DBL_EPSILON) {
    vec_sub(x[frag->atoms[iatom]], x[frame->ref_atom1], a);
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

// Pre-calculate modal, directional DoF of rigid fragments, since these can just be rotated
// with the fragment rather than re-calculating each time
void dofulator_precalculate_rigid(Dofulator ctx, const double* mass, const double (*x)[3]) {
  assert(ctx);
  ctx->rigid_ref_frames = realloc(ctx->rigid_ref_frames, sizeof(RefFrame) * ctx->rigid_frags.n_fragments);
  // TODO: better error handling
  assert(ctx->rigid_ref_frames);
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
        rij[0] = x[atom][0] - x[root][0];
        rij[1] = x[atom][1] - x[root][1];
        rij[2] = x[atom][2] - x[root][2];

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
      fragment_solve_dof(ctx, frag, mJ, Q);
      fragment_set_frame(frag, &ctx->rigid_ref_frames[ifrag], x);
    }
  }
}

// Calculate rotation quaternion for each rigid fragment.
void dofulator_calculate_rigid(Dofulator ctx, const double (*x)[3]) {
  for (size_t ifrag = 0; ifrag < ctx->rigid_frags.n_fragments; ++ifrag) {
    ctx->rigid_ref_frames[ifrag].current_rot = frame_get_rotation(&ctx->rigid_ref_frames[ifrag], x);
  }
}

// Calculate dof_total and dof for all semirigid fragments based on given masses
// and atom locations (indexed by atom index)
void dofulator_calculate_semirigid(Dofulator ctx, const double* mass, const double (*x)[3]) {
  for (
    size_t ifrag = 0, *n_frags = ctx->n_semirigid;
    n_frags < ctx->n_semirigid + ctx->max_semirigid_atoms;
    ++n_frags
  ) {
    if (*n_frags == 0) { continue; }
    const size_t n_atoms = ctx->semirigid_frags.fragments[ifrag].n_atoms;
    const size_t mat_size = 3*n_atoms * 3*n_atoms;
    double* mJ = ctx->working_buf;
    double* Q = mJ + mat_size;
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
          memcpy(&((*J)[3*i][0]), &((*J)[3*pred_i][0]), sizeof(double) * 3*n_atoms * 3);

          // Get vector from predecessor to atom i
          rij[0] = x[atom][0] - x[pred][0];
          rij[1] = x[atom][1] - x[pred][1];
          rij[2] = x[atom][2] - x[pred][2];

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
        double (*K)[frag->n_loops][3*n_atoms] = (double(*)[frag->n_loops][n_modes])mJ;
        for (size_t k = 0; k < frag->n_loops; ++k) {
          const Bond b = frag->loop_closures[k];
          const size_t i = ctx->atom_frag_idx[b.ai];
          const size_t j = ctx->atom_frag_idx[b.aj];

          // T^T is a row vector for bonds
          rij[0] = x[b.aj][0] - x[b.ai][0];
          rij[1] = x[b.aj][1] - x[b.ai][1];
          rij[2] = x[b.aj][2] - x[b.ai][2];

          // Calculate J_j - J_i, store 3 x 3*n_atoms result in Q
          memcpy(Q, &((*J)[3*j][0]), sizeof(double) * n_modes*3);
          cblas_daxpy(3*n_modes, -1.0, &((*J)[3*i][0]), 1, Q, 1);
          
          // Calculate row k of K = T^T (J_j - J_i)
          // Use (J_j - J_i)^T x T and store result in K[k][:]
          cblas_dgemv(CblasRowMajor, CblasTrans, 3, n_modes, 1.0, Q, n_modes, rij, 1, 0., &((*K)[k][0]), 1);
        }

        // Find null(K) from low singular values of K
        double work_query;
        double* S = frag->dof_total;  // Temporarily store singular values in dof_total
        lapack_int info = LAPACKE_dgesvd_work(
          LAPACK_ROW_MAJOR,
          'N', // Don't calculate any of U
          'A', // Calculate all of V^T
          n_loops, n_modes,
          &((*K)[0][0]), n_loops,
          S,
          NULL, frag->n_loops,
          VT, n_modes,
          &work_query,
          -1
        );
        // TODO: better error handling
        assert(info == 0);
        // Need more working space, so allocate it
        if (work_query > ctx->svd_working_size) {
          ctx->svd_working_size = (lapack_int)work_query;
          ctx->svd_working = realloc(ctx->svd_working, sizeof(double) * ctx->svd_working_size);
          // TODO: better error handling
          assert(ctx->svd_working);
        }
        info = LAPACKE_dgesvd_work(
          LAPACK_ROW_MAJOR,
          'N', // Don't calculate any of U
          'A', // Calculate all of V^T
          n_loops, n_modes,
          &((*K)[0][0]), n_loops,
          S,
          NULL, frag->n_loops,
          VT, n_modes,
          ctx->svd_working,
          ctx->svd_working_size
        );
        // TODO: better error handling
        assert(info == 0);

        // Find rank of null(K)
        double thresh = ctx->null_space_thresh * frag->dof_total[0];   // Singular values in descending order
        size_t rank_K = 0;
        while (S[rank_K] > thresh) { ++rank_K; }

        // Project Jacobian onto null(K).
        // null(K) = VT[modes:3*n_atoms][..]
        cblas_dgemm(
          CblasRowMajor, CblasNoTrans, CblasTrans,
          n_modes, n_modes - rank_K, n_modes,
          1.0, &((*J)[0][0]), n_modes,
          VT, n_modes,
          0.0, mJ, n_modes
        );

        frag->n_modes = n_modes - rank_K;
      }

      // Scale rows by sqrt(mass) of corresponding atom
      for (AtomTag a = 0; a < n_atoms; ++a) {
        // TODO: pre-calculate sqrt(mass) and store on fragment?
        double root_mass = sqrt(mass[frag->atoms[a]]);
        cblas_dscal(frag->n_modes, root_mass, &mJ[3*a*3*n_atoms], 1);
        cblas_dscal(frag->n_modes, root_mass, &mJ[(3*a+1)*3*n_atoms], 1);
        cblas_dscal(frag->n_modes, root_mass, &mJ[(3*a+2)*3*n_atoms], 1);
      }

      fragment_solve_dof(ctx, frag, mJ, Q);
    }
  }
}

// Calculate DoF of all fragments in the context
void dofulator_calculate(Dofulator ctx, const double* mass, const double (*x)[3]) {
  dofulator_calculate_rigid(ctx, x);
  dofulator_calculate_semirigid(ctx, mass, x);
}

// Get the directional DoF of the atom with index `atom_idx`
void dofulator_get_dof_atom_directional(const struct Dofulator* ctx, AtomTag atom_idx, double dof[3]) {
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
    quat_rotate_vec(ctx->rigid_ref_frames[frag_idx.idx].current_rot, x);
    quat_rotate_vec(ctx->rigid_ref_frames[frag_idx.idx].current_rot, y);
    quat_rotate_vec(ctx->rigid_ref_frames[frag_idx.idx].current_rot, z);
    dof[0] = frag->dof_total[3*i]*x[0]*x[0] + frag->dof_total[3*i+1]*y[0]*y[0] + frag->dof_total[3*i+2]*z[0]*z[0];
    dof[1] = frag->dof_total[3*i]*x[1]*x[1] + frag->dof_total[3*i+1]*y[1]*y[1] + frag->dof_total[3*i+2]*z[1]*z[1];
    dof[2] = frag->dof_total[3*i]*x[2]*x[2] + frag->dof_total[3*i+1]*y[2]*y[2] + frag->dof_total[3*i+2]*z[2]*z[2];
  } else {
    dof[0] = frag->dof_total[3*i];
    dof[1] = frag->dof_total[3*i+1];
    dof[2] = frag->dof_total[3*i+2];
  }

}

// Get the total DoF of the atom with index `atom_idx`
double dofulator_get_dof_atom(const struct Dofulator* ctx, AtomTag atom_idx) {
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
