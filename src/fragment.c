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
 * Allocate space for a fragment and set pointers
 *
*/
static inline Fragment_t* fragment_alloc(AtomList atoms) {
  if (atoms.n == 0) return NULL;

  // Need to store three matrices each up to 3n x 3n,
  // plus one array of length n
  unsigned mat_sz = 3 * atoms.n * 3 * atoms.n;
  Fragment_t* frag = (Fragment_t*)calloc(1,
    sizeof(Fragment_t) + sizeof(double) * (3*mat_sz + atoms.n));
  if (!frag) return NULL;

  frag->natoms = atoms.n;
  frag->mC = &frag->R[mat_sz];
  frag->mCR = &frag->R[2*mat_sz];
  frag->Imodal = &frag->R[3*mat_sz];

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

  double (*mC)[3*atoms.n][3*atoms.n] = (double(*)[3*atoms.n][3*atoms.n])frag->mC;
  Vec3 rij;
  (*mC)[0][0] = (*mC)[1][1] = (*mC)[2][2] = sqrt(atoms.mass[0]);
  for (unsigned i = 1; i < atoms.n; ++i) {
    const double root_m = sqrt(atoms.mass[i]);
    (*mC)[3*i][0] = (*mC)[3*i+1][1] = (*mC)[3*i+2][2] = root_m;
    rij = vec_sub(&atoms.pos[3*i], &atoms.pos[0]);

    /* Store negative cross operator matrix
     *  0   z   -y
     * -z   0    x
     *  y  -x    0
     */
    (*mC)[3*i][4] = root_m * rij.z;
    (*mC)[3*i][5] = -root_m * rij.y;

    (*mC)[3*i+1][3] = -root_m * rij.z;
    (*mC)[3*i+1][5] = root_m * rij.x;

    (*mC)[3*i+2][3] = root_m * rij.y;
    (*mC)[3*i+2][4] = -root_m * rij.x;
  }

  // Maximum possible number of modes
  frag->nmodes = atoms.n > 1 ? 6 : 3;

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
  const unsigned n = frag->natoms;
  const unsigned modes = frag->nmodes;

  // Get lab-frame inertia tensor == C^T C
  // Store in R since it will be overwritten by eigenvectors
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, modes, modes, 3*n,
              1.0, frag->mC, 3*n, frag->mC, 3*n, 0.0, frag->R, modes);

  // Find eigenvalues/eigenvectors of inertia tensor
  lapack_int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'L', modes,
                                  frag->R, modes, frag->Imodal);
  if (info) {
    fragment_destroy((Fragment)frag);
    return NULL;
  }

  // Calculate modal inertias:
  // I_mode = R^T C^T M C R
  //        = R_mode^T lambda_mode R_mode
  //        = lambda_mode R_mode dot R_mode
  for (unsigned m = 0; m < modes; ++m) {
    frag->Imodal[m] *= cblas_ddot(modes, &frag->R[m], modes, &frag->R[m], modes);
  }

  // Store the full transform M^1/2 C R for later to get per-atom modal inertia
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3*n, modes, modes,
              1.0, frag->mC, 3*n, frag->R, modes, 0.0, frag->mCR, modes);

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
  const unsigned n = frag->natoms;
  const unsigned modes = frag->nmodes;
  double (*mCR)[3*n][modes] = (double(*)[3*n][modes])frag->mCR;

  // DoF == sum_modes <atom modal inertia> / <modal inertia>
  // Atom modal directional inertias == diag(mCR^T mCR) == diag(R^T C^T M C R)

  for (unsigned m = 0; m < modes; ++m) {
    if (frag->Imodal[m] > 100*DBL_EPSILON) {
      double Ix = (*mCR)[3*atom][m] * (*mCR)[3*atom][m];
      double Iy = (*mCR)[3*atom+1][m] * (*mCR)[3*atom+1][m];
      double Iz = (*mCR)[3*atom+2][m] * (*mCR)[3*atom+2][m];
      dof += (Ix + Iy + Iz) / frag->Imodal[m];
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
  const unsigned n = frag->natoms;
  const unsigned modes = frag->nmodes;
  double (*mCR)[3*n][modes] = (double(*)[3*n][modes])frag->mCR;

  // DoF == sum_modes <atom modal inertia> / <modal inertia>
  // Atom modal directional inertias == diag(mCR^T mCR) == diag(R^T C^T M C R)
  // Dot product with direction vector to get inertia in that direction.

  for (unsigned m = 0; m < modes; ++m) {
    if (frag->Imodal[m] > 100*DBL_EPSILON) {
      double Ix = (*mCR)[3*atom][m] * (*mCR)[3*atom][m];
      double Iy = (*mCR)[3*atom+1][m] * (*mCR)[3*atom+1][m];
      double Iz = (*mCR)[3*atom+2][m] * (*mCR)[3*atom+2][m];
      dof += ( Ix*dir[0] + Iy*dir[1] + Iz*dir[2] ) / frag->Imodal[m];
    }
  }
  return dof;
}
