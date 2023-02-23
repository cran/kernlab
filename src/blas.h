#ifndef KERNLAB_BLAS_H
#define KERNLAB_BLAS_H

// For forward compatibility with >= R 4.2.0
#ifndef USE_FC_LEN_T
# define USE_FC_LEN_T
#endif

#include <R_ext/BLAS.h>

// For backwards compatibility with <= R 3.6.2
// where `FCONE` wasn't defined by R yet
#ifndef FCONE
# define FCONE
#endif

#endif
