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
static inline Fragment_t* fragment_alloc(AtomList atoms, unsigned max_modes) {
  if (atoms.n == 0) return NULL;

  // Need to store two matrices up to 3n x max_modes,
  // one at max_modes x max_modes,
  // plus one array of length max_modes
  unsigned mat_sz = 3 * atoms.n * max_modes;
  unsigned mode_mat_sz = max_modes * max_modes;
  Fragment_t* frag = (Fragment_t*)calloc(1,
    sizeof(Fragment_t) + sizeof(double) * (mode_mat_sz + 2*mat_sz + max_modes));
  if (!frag) return NULL;

  frag->natoms = atoms.n;
  frag->mC     = &frag->R[mode_mat_sz];
  frag->mCR    = &frag->R[mode_mat_sz + mat_sz];
  frag->Imodal = &frag->R[mode_mat_sz + 2*mat_sz];

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
 * Bond list will be sanitized (mutated)
 *
 * WARNING: Currently cannot handle rings - only tree topology!
 *
*/
Fragment fragment_create_semirigid(AtomList atoms, BondList bonds) {
  assert(bonds.n > 0);

  // Make sure there are no loops, assuming all atoms connected
  assert(atoms.n - bonds.n == 1);

  // Make sure bonds are sorted, and bond.i < bond.j
  sanitize_bonds(bonds);

  // Maximum possible number of modes
  unsigned max_modes = 3*atoms.n;
  Fragment_t* frag = fragment_alloc(atoms, max_modes);
  if (!frag) return NULL;

  frag->nmodes = max_modes;

  double (*mC)[3*atoms.n][max_modes] = (double(*)[3*atoms.n][max_modes])frag->mC;
  Vec3 rij;

  // Add atoms breadth-first so that rows can be copied as new atoms are encountered
  unsigned base = bonds.bonds[0].i;
  (*mC)[3*base][3*base] = 1.0;
  (*mC)[3*base+1][3*base+1] = 1.0;
  (*mC)[3*base+2][3*base+2] = 1.0;
  for (unsigned b = 0; b < bonds.n; ++b) {
    const unsigned i = bonds.bonds[b].i;
    const unsigned j = bonds.bonds[b].j;
    rij = vec_sub(&atoms.pos[3*j], &atoms.pos[3*i]);

    // Copy rows from i to get previous connectivity
    memcpy((*mC)[3*j], (*mC)[3*i], sizeof(double)*3*max_modes);

    /* Store negative cross operator matrix
     *  0   z   -y
     * -z   0    x
     *  y  -x    0
     */
    (*mC)[3*j][3*j+1] = rij.z;
    (*mC)[3*j][3*j+2] = -rij.y;

    (*mC)[3*j+1][3*j] = -rij.z;
    (*mC)[3*j+1][3*j+2] = rij.x;

    (*mC)[3*j+2][3*j] = rij.y;
    (*mC)[3*j+2][3*j+1] = -rij.x;
  }

  // Apply mass scaling
  for (unsigned i = 0; i < atoms.n; ++i) {
    const double root_m = sqrt(atoms.mass[i]);
    for (unsigned d = 0; d < 3; ++d) {
      for (unsigned m = 0; m < max_modes; ++m) {
       (*mC)[3*i+d][m] *= root_m;
      }
    }
  }

  return (Fragment)fragment_eval(frag);
}


/* =============================================================================
 *
 * Create a rigid fragment from an atom list.
 * All atoms treated as part of a single rigid body.
 *
*/
Fragment fragment_create_rigid(AtomList atoms) {
  // Maximum possible number of modes
  unsigned max_modes = atoms.n > 1 ? 6 : 3;
  Fragment_t* frag = fragment_alloc(atoms, max_modes);
  if (!frag) return NULL;

  frag->nmodes = max_modes;

  double (*mC)[3*atoms.n][max_modes] = (double(*)[3*atoms.n][max_modes])frag->mC;
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


  return (Fragment)fragment_eval(frag);
}


/* =============================================================================
 *
 * Evaluate fragment properties.
 * Assumes mC[3*natoms][nmodes] has been populated.
 * R and mCR will be calculated by this routine.
 *
*/
Fragment_t* fragment_eval(Fragment_t* frag) {
  assert(frag);

  const unsigned n = frag->natoms;
  const unsigned modes = frag->nmodes;

  // Get lab-frame inertia tensor == C^T C
  // Store in R since it will be overwritten by eigenvectors
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, modes, modes, 3*n,
              1.0, frag->mC, modes, frag->mC, modes, 0.0, frag->R, modes);

  // Find eigenvalues/eigenvectors of inertia tensor
  lapack_int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'L', modes,
                                  frag->R, modes, frag->Imodal);
  if (info) {
    fragment_destroy((Fragment)frag);
    return NULL;
  }

  // Modal inertia = R^T C^T M C R
  // For mode m:
  //           I_m = lambda_m R_m^T R_m
  //               = lambda_m

  // Store the full transform M^1/2 C R for later to get per-atom modal inertia
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3*n, modes, modes,
              1.0, frag->mC, modes, frag->R, modes, 0.0, frag->mCR, modes);

  return frag;
}

/* =============================================================================
 *
 * Calculate the total DoF of the fragment
 *
*/
double fragment_dof(const Fragment fragment) {
  const Fragment_t* frag = (const Fragment_t*)fragment;
  double dof = 0.0;
  for (unsigned i = 0; i < frag->natoms; ++i) {
    dof += fragment_dof_atom(fragment, i);
  }
  return dof;
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
double fragment_dof_atom_dir(Fragment fragment, unsigned atom, const double dir[3]) {
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


/* =============================================================================
 *
 * Comparison routine for sorting bond list.
 *
*/
static int cmp_bonds(const void* lhs, const void* rhs) {
  Bond* l = (Bond*)lhs;
  Bond* r = (Bond*)rhs;
  int di = l->i - r->i;
  if (di) return l->j - r->j;
  return di;
}


/* =============================================================================
 *
 * Sanitize the bond list so that i < j and list is sorted by i, then j
 *
*/
static inline void sanitize_bonds(BondList bonds) {
  // Swap so that i < j
  for (unsigned b = 0; b < bonds.n; ++b) {
    assert(bonds.bonds[b].i != bonds.bonds[b].j);
    if (bonds.bonds[b].j < bonds.bonds[b].i) {
      unsigned tmp = bonds.bonds[b].i;
      bonds.bonds[b].i = bonds.bonds[b].j;
      bonds.bonds[b].j = tmp;
    }
  }

  qsort(bonds.bonds, bonds.n, sizeof(Bond), cmp_bonds);

  for (unsigned b = 1; b < bonds.n; ++b) {
    assert(memcmp(&bonds.bonds[b], &bonds.bonds[b-1], sizeof(Bond)));
  }
}
