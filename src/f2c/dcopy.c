#include "blas.h"
#include <string.h>  /* Needed for memcpy() */

int dcopy_(long *n, double *sx, long *incx, double *sy, long *incy)
{
  long i, ix, iy, nn, iincx, iincy;

  /* copies a vector, x, to a vector, y.   
     uses unrolled loops for increments equal to 1.   
     jack dongarra, linpack, 3/11/78.   
     modified 12/3/93, array(1) declarations changed to array(*) */

  /* Dereference inputs */
  nn = *n;
  iincx = *incx;
  iincy = *incy;

  if (nn > 0)
  {
    if (iincx == 1 && iincy == 1) /* code for both increments equal to 1 */
    {
      memcpy( sy, sx, nn * sizeof(*sy) );
    }
    else /* code for unequal increments or equal increments not equal to 1 */
    {
      ix = iincx >= 0 ? 0 : (1 - nn) * iincx;
      iy = iincy >= 0 ? 0 : (1 - nn) * iincy;
      for (i = 0; i < nn; i++)
      {
        sy[iy] = sx[ix];
        ix += iincx;
        iy += iincy;
      }
    }
  }

  return 0;
} /* dcopy_ */
