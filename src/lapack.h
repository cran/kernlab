#ifndef KERNLAB_LAPACK_H
#define KERNLAB_LAPACK_H

// For forward compatibility with >= R 4.2.0
#ifndef USE_FC_LEN_T
# define USE_FC_LEN_T
#endif

#include <R_ext/Lapack.h>

// For backwards compatibility with <= R 3.6.2
// where `FCONE` wasn't defined by R yet
#ifndef FCONE
# define FCONE
#endif

#endif
