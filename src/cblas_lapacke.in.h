#include <@INCLUDE_CBLAS@>
#include <@INCLUDE_LAPACKE@>

#ifdef BLA_VENDOR_OPENBLAS

// OpenBLAS provides dsum as an extension, so use that if available
#define dsum(N, in) cblas_dsum((N), (in), 1)

#else

#include <stddef.h>

/*******************************************************************************
 * Portable dsum impl which encourages auto-vectorisation without -ffast-math.
 * Better than naive loop, but could be improved with non-portable intrinsics.
*/
static inline double dsum(size_t N, const double* in) {
  double buf[8] = {0.};
  size_t m = 0, i;

  for (; m + 8 < N; m += 8) {
    for (i = 0; i < 8; ++i) {
      buf[i] += in[m + i];
    }
  }

  for (i = 0; i < N - m; ++i) {
    buf[i] += in[m + i];
  }

  buf[0] += buf[4];
  buf[1] += buf[5];
  buf[2] += buf[6];
  buf[3] += buf[7];

  buf[0] += buf[2];
  buf[1] += buf[3];

  return buf[0] + buf[1];
}

#endif
