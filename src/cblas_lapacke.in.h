#include <@INCLUDE_CBLAS@>
#include <@INCLUDE_LAPACKE@>

#ifdef SCIPY_OPENBLAS
#define cblas(fn) scipy_cblas_##fn
#else
#define cblas(fn) cblas_##fn
#endif
