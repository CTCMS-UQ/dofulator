#ifdef _MSC_VER
// msvc doesn't support function arguments declared as foo[restrict]
// Use macro to preserve size information in declaration
#define NOALIAS_ARR(name, size) *restrict name
#else
#define NOALIAS_ARR(name, size) name[restrict size]
#endif
