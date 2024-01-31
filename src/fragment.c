#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <openblas/cblas.h>
#include <openblas/lapacke.h>
#include <string.h>

#include "fragment.h"
#include "vec3.h"

/* =============================================================================
 *
 * Allocate space for a fragment and initialise its atom list to `atoms` and 
 * its mass matrix (`M`) from the atom masses.
 *
*/
static inline Fragment_t* fragment_alloc(AtomList atoms) {
  if (atoms.n == 0) return NULL;

  // Need to store five matrices each up to 3n x 3n,
  // plus one more as scratch space
  unsigned mat_sz = 3 * atoms.n * 3 * atoms.n;
  Fragment_t* frag = (Fragment_t*)calloc(1,
    sizeof(Fragment_t) + sizeof(double) * 6*mat_sz);
  if (!frag) return NULL;

  frag->atoms = atoms;
  // buf[0..mat_sz] used as temporary scratch space
  frag->M = &frag->buf[mat_sz];
  frag->C = &frag->buf[2*mat_sz];
  frag->R = &frag->buf[3*mat_sz];
  frag->Imodal  = &frag->buf[4*mat_sz];
  frag->Iatom = &frag->buf[5*mat_sz];

  // Construct mass matrix - block diagonal of atomic mass times 3x3 identity
  double (*M)[3*atoms.n][3*atoms.n] = (double(*)[3*atoms.n][3*atoms.n])frag->M;
  for (unsigned i = 0; i < atoms.n; ++i) {
    (*M)[3*i][3*i]     = atoms.mass[i];
    (*M)[3*i+1][3*i+1] = atoms.mass[i];
    (*M)[3*i+2][3*i+2] = atoms.mass[i];
  }

  return frag;
}


/* =============================================================================
 *
 * Clean up a fragment.
 *
*/
void fragment_destroy(Fragment fragment) {
  free(fragment);
}


/* =============================================================================
 *
 * Create a semi-rigid fragment from an atom list and bond list.
 * Treats all bonds as fixed length with flexible angles, and assumes all atoms
 * are connected in a single graph.
 *
 * WARNING: Currently cannot handle rings - only tree topology!
 *
*/
Fragment fragment_create_semirigid(AtomList atoms, BondList bonds) {
  Fragment_t* frag = fragment_alloc(atoms);
  if (!frag) return NULL;

  return (Fragment)frag;
}


/* =============================================================================
 *
 * Create a rigid fragment from an atom list.
 * All atoms treated as part of a single rigid body.
 *
*/
Fragment fragment_create_rigid(AtomList atoms) {
  Fragment_t* frag = fragment_alloc(atoms);
  if (!frag) return NULL;

  double (*C)[3*atoms.n][3*atoms.n] = (double(*)[3*atoms.n][3*atoms.n])frag->C;
  Vec3 rij;
  (*C)[0][0] = (*C)[1][1] = (*C)[2][2] = 1.0;
  for (unsigned i = 1; i < atoms.n; ++i) {
    (*C)[3*i][0] = (*C)[3*i+1][1] = (*C)[3*i+2][2] = 1.0;
    rij = vec_sub(&atoms.pos[3*i], &atoms.pos[0]);

    /* Store negative cross operator matrix
     *  0   z   -y
     * -z   0    x
     *  y  -x    0
     */
    (*C)[3*i][4] = rij.z;
    (*C)[3*i][5] = -rij.y;

    (*C)[3*i+1][3] = -rij.z;
    (*C)[3*i+1][5] = rij.x;

    (*C)[3*i+2][3] = rij.y;
    (*C)[3*i+2][4] = -rij.x;
  }

  // Maximum possible number of modes
  frag->n_modes = atoms.n > 1 ? 6 : 3;

  return (Fragment)fragment_eval(frag);
}


/* =============================================================================
 *
 * Evaluate fragment properties.
 * Assumes C[n_modes][3*atoms.n] and M[3*atoms.n][3*atoms.n] are populated.
 * All others will be calculated by this routine.
 *
*/
Fragment_t* fragment_eval(Fragment_t* frag) {
  const unsigned n = frag->atoms.n;
  const blasint modes = frag->n_modes;

  // Get lab-frame inertia tensor
  // TODO: better blas routines for this? Take advantage of M being diagonal?
  // C^T M
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, modes, 3*n, 3*n,
              1.0, frag->C, 3*n, frag->M, 3*n, 0.0, frag->buf, 3*n);
  // (C^T M) C
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, modes, modes, 3*n,
              1.0, frag->buf, 3*n, frag->C, 3*n, 0.0, frag->R, modes);
  // Stored total lab-frame inertia in R since will be overwritten with eigenvectors
  // and won't be needed after that point.

  // Find eigenvalues/eigenvectors of inertia tensor
  lapack_int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'L', modes,
                                  frag->R, modes, frag->Imodal);
  if (info) {
    fragment_destroy((Fragment)frag);
    return NULL;
  }

  // Calculate modal inertias:
  // I_mode = R^T C^T M C R = R_mode^T lambda_mode R_mode = lambda_mode R_mode dot R_mode
  for (unsigned m = 0; m < frag->n_modes; ++m) {
    frag->Imodal[m] *= cblas_ddot(modes, &frag->R[m], modes, &frag->R[m], modes);
  }

  // The transform from modal motion to lab motion == C R
  // No longer need buf, so store R <- C R
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3*n, modes, modes,
              1.0, frag->C, 3*n, frag->R, modes, 0.0, frag->buf, modes);

  // For each atom, store the x,y,z modal inertia in Iatom
  double (*Iatom)[modes][3*n] = (double(*)[modes][3*n])frag->Iatom;
  double (*CR)[3*n][modes] = (double(*)[3*n][modes])frag->buf;

  for (unsigned i = 0; i < n; ++i) {
    for (unsigned m = 0; m < frag->n_modes; ++m) {
      (*Iatom)[m][3*i]   = frag->atoms.mass[i] * (*CR)[3*i][m]   * (*CR)[3*i][m];
      (*Iatom)[m][3*i+1] = frag->atoms.mass[i] * (*CR)[3*i+1][m] * (*CR)[3*i+1][m];
      (*Iatom)[m][3*i+2] = frag->atoms.mass[i] * (*CR)[3*i+2][m] * (*CR)[3*i+2][m];
    }
  }

  return frag;
}


/* =============================================================================
 *
 * Calculate the DoF of atom with index `atom`
 *
*/
double fragment_dof_atom(const Fragment fragment, unsigned atom) {
  const Fragment_t* frag = (const Fragment_t*)fragment;
  double dof = 0.0;
  const unsigned n = frag->atoms.n;
  const unsigned modes = frag->n_modes;
  double (*Iatom)[modes][3*n] = (double(*)[modes][3*n])frag->Iatom;
  for (unsigned m = 0; m < modes; ++m) {
    if (frag->Imodal[m] > 100*DBL_EPSILON) {
      dof += ( (*Iatom)[m][3*atom] + (*Iatom)[m][3*atom+1]
             + (*Iatom)[m][3*atom+2] ) / frag->Imodal[m];
    }
  }
  return dof;
}

/* =============================================================================
 *
 * Calculate the DoF of atom with index `atom` in direction `dir`.
 * `dir` must be a unit vector.
 *
*/
double fragment_dof_atom_dir(Fragment fragment, unsigned atom, double dir[3]) {
  assert(fabs(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2] - 1.0) < 100*DBL_EPSILON);
  const Fragment_t* frag = (const Fragment_t*)fragment;
  double dof = 0.0;
  const unsigned n = frag->atoms.n;
  const unsigned modes = frag->n_modes;
  double (*Iatom)[modes][3*n] = (double(*)[modes][3*n])frag->Iatom;
  for (unsigned m = 0; m < modes; ++m) {
    if (frag->Imodal[m] > 100*DBL_EPSILON) {
      dof += ( (*Iatom)[m][3*atom]   * dir[0]
             + (*Iatom)[m][3*atom+1] * dir[1]
             + (*Iatom)[m][3*atom+2] * dir[2]
             ) / frag->Imodal[m];
    }
  }
  return dof;
}
