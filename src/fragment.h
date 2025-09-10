#ifndef DOFULATOR_FRAGMENT_H
#define DOFULATOR_FRAGMENT_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "dofulator.h"
#include "quaternion.h"

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

/*******************************************************************************
 * Single fragment
*/
struct Fragment {
  uint32_t n_atoms;         // Number of atoms in the fragment
  uint32_t n_modes;         // Maximum number of modes this fragment could have
  uint32_t rigid : 1;       // 1 = rigid, 0 = semirigid
  uint32_t invalid : 1;     // 1 = fragment has been invalidated, in which case root_atom
                            // is the index of the fragment it was merged into
  uint32_t n_loops : 30;    // Number of kinematic loops (length of `loop_closures`)

  AtomTag root_atom;        // Root atom of the fragment's tree.
                            // For merged fragments, this is the idx of the parent
                            // fragment in the context's `frag_map`.

  AtomTag* atoms;           // Indices of atoms in this fragment. Stored in ascending order.
  Bond* loop_closures;      // Rigid bonds that create kinematic loops

  double* dof;              // 3*n_atoms x n_modes matrix of M^1/2 J Q
  double* dof_total;        // 3*n_atoms vector of total direciontal dof. Rotational only for rigid bodies.
  double* dof_trans;        // n_atoms elements for translational DoF of each atom (mass frac.)
};


/*******************************************************************************
 * Managed list of fragments
*/
typedef struct FragmentList {
  size_t capacity;      // Current capacity
  size_t n_fragments;   // Number of fragments
  Fragment* fragments;  // Individual fragments
} FragmentList;


/*******************************************************************************
 * Reference frame for rigid fragments
*/
typedef struct RefFrame {
  size_t ref_atom1, ref_atom2, ref_atom3; // Atoms to use to construct the reference frame.
                                          // ref_atom2 == ref_atom3 if only 2 atoms in frag
  double r12[3], r13_perp[3];   // r12 and (r12 x r13) vectors in the original coordinates
  Quaternion current_rot;       // Current rotation relative to the original frame.
} RefFrame;

/*******************************************************************************
 * Mapping for atom onto fragment
*/
typedef struct FragIndex {
  size_t has_frag : 1;                  // 0 = free atom, 1 = atom is in a fragment
  size_t rigid : 1;                     // 0 = semi-rigid fragment, 1 = rigid
  size_t idx : (sizeof(size_t)*8 - 2);  // Index of fragment in the relevant list
} FragIndex;


/*******************************************************************************
 * Create a rigid or semi-rigid fragment with the specified number of atoms
*/
static inline Fragment fragment_create(AtomTag n_atoms, bool rigid) {
  return (Fragment){
    .n_atoms = n_atoms,
    .rigid = rigid,
    .n_loops = 0,
  };
}


/*******************************************************************************
 * Clear a FragmentList and free its buffer
*/
static void fragmentlist_destroy(FragmentList* self) {
  if (!self) return;
  if (self->fragments) {
    for (Fragment* frag = self->fragments; frag < self->fragments + self->n_fragments; ++frag) {
      if (frag->loop_closures) { free(frag->loop_closures); }
    }
    free(self->fragments);
  }
  self->n_fragments = 0;
  self->capacity = 0;
}


/*******************************************************************************
 * Set the maximum number of expected modes of a finished fragment.
 * 6 for rigid bodies, 3 * n_atoms for semi-rigid.
*/
static inline void fragment_set_max_modes(Fragment* frag) {
  frag->n_modes = frag->rigid ? 3 : 3 * frag->n_atoms;
}


/*******************************************************************************
 * Add a new rigid fragment with `n_atoms` atoms to `list`.
*/
#define fragmentlist_add_rigid(list, n_atoms) fragmentlist_add_new((list), (n_atoms), true)

/*******************************************************************************
 * Add a new semirigid fragment with `n_atoms` atoms to `list`.
*/
#define fragmentlist_add_semirigid(list, n_atoms) fragmentlist_add_new((list), (n_atoms), false)

/*******************************************************************************
 * Add a new fragment with `n_atoms` atoms to `list`.
 * `rigid` should be `true` for rigid fragments, or `false` for semirigid.
*/
typedef struct IndexResult { DofulatorResult status; size_t idx; } IndexResult;
static IndexResult fragmentlist_add_new(FragmentList* list, size_t n_atoms, bool rigid) {
  // Find the next empty slot
  // Allocate if needed
  if (list->n_fragments >= list->capacity) {
    // List needs to grow, so pick a new capacity and extend it
    list->capacity = list->capacity < 4 ? 4 : list->capacity * 2;
    Fragment* new = realloc(list->fragments, sizeof(Fragment) * list->capacity);
    if (unlikely(!new)) { return (IndexResult){DOF_ALLOC_FAILURE, 0}; }
    list->fragments = new;
  }

  // Create the new fragment.
  // Allocation and filling of atom list/buffers comes later.
  size_t idx = list->n_fragments++;
  list->fragments[idx] = fragment_create(n_atoms, rigid);

  return (IndexResult){DOF_SUCCESS, idx};
}


/*******************************************************************************
 * Get index of the top-level fragment containing fragment `idx`.
 * Used to get updated fragment index from an invalidated (merged) fragment.
*/
static size_t fragmentlist_get_idx(FragmentList fraglist, size_t idx) {
  if (!fraglist.fragments[idx].invalid) {
    return idx;
  }
  return fragmentlist_get_idx(fraglist, fraglist.fragments[idx].root_atom);
}


/*******************************************************************************
 * Get the index of the fragment which contains the atom with index `atom_idx`.
 * May mutate `frag_map` if it points to an invalidated (merged) fragment.
 * Assumes fraglist is the correct one to be indexed by `frag_map` for the given atom.
*/
static size_t fragmentlist_get_fragment_idx(FragmentList fraglist, FragIndex* frag_map, AtomTag atom_idx) {
  size_t guess = frag_map[atom_idx].idx;
  size_t idx = fragmentlist_get_idx(fraglist, guess);
  if (idx != guess) {
    frag_map[atom_idx].idx = idx;
  }
  return idx;
}


#endif
