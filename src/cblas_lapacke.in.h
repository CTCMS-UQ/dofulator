#include <@INCLUDE_CBLAS@>
#include <@INCLUDE_LAPACKE@>

#ifdef SCIPY_OPENBLAS
#define cblas(fn) scipy_cblas_##fn
#define lapack(fn) scipy_##fn##_
#else
#define cblas(fn) cblas_##fn
#define lapack(fn) fn##_
#endif

// Wrappers for fortran LAPACK calls
#ifdef LAPACK_FORTRAN_STRLEN_END
#define LAPACK_FN_END , FORTRAN_STRLEN, FORTRAN_STRLEN
#define LAPACK_CALL_END , 1, 1
#else
#define LAPACK_FN_END
#define LAPACK_CALL_END
#endif

extern void lapack(dsyev)(
  const char* jobz, const char* uplo, const lapack_int* n, double* a, const lapack_int* lda,
  double* w, double* work, const lapack_int* lwork, lapack_int* info
  LAPACK_FN_END
);
static lapack_int fLAPACK_dsyev(char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w, double* work, lapack_int lwork) {
  lapack_int info = 0;
  lapack(dsyev)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info LAPACK_CALL_END);
  return info < 0 ? info - 1 : info;
}

extern void lapack(dgesvd)(
  const char* jobu, const char* jobvt, const lapack_int* m, const lapack_int* n,
  double* a, const lapack_int* lda, double* s, double* u, const lapack_int* ldu,
  double* vt, const lapack_int* ldvt, double* work, const lapack_int* lwork, lapack_int* info
  LAPACK_FN_END
);
static lapack_int fLAPACK_dgesvd(char jobu, char jobvt, lapack_int m, lapack_int n, double* a, lapack_int lda, double* s, double* u, lapack_int ldu, double* vt, lapack_int ldvt, double* work, lapack_int lwork) {
  lapack_int info = 0;
  lapack(dgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info LAPACK_CALL_END);
  return info < 0 ? info - 1 : info;
}
