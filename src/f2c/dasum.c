#include "blas.h"
#include <math.h> /* needed for fabs() */

double dasum_(long *n, double *sx, long *incx)
{
  long i, m, nn, iincx;
  double stemp;

/* takes the sum of the absolute values.   
   uses unrolled loops for increment equal to one.   
   jack dongarra, linpack, 3/11/78.   
   modified 3/93 to return if incx .le. 0.   
   modified 12/3/93, array(1) declarations changed to array(*) */

  /* Dereference inputs */
  nn = *n;
  iincx = *incx;

  stemp = 0.0;
  if (nn > 0 && iincx > 0)
  {
    if (iincx == 1) /* code for increment equal to 1 */
    {
      m = nn % 6;
      for (i = m; i < nn; i += 6)
      {
        stemp += fabs(sx[i]) + fabs(sx[i+1]) + fabs(sx[i+2]) +
          fabs(sx[i+3]) + fabs(sx[i+4]) + fabs(sx[i+5]);
      }
      if (m != 0) /* clean-up loop */
      {
        for (i = 0; i < m; ++i)
        {
          stemp += fabs(sx[i]);
        }
      }
    }
    else /* code for increment not equal to 1 */
    {
      for (i=(nn-1)*iincx; i>=0; i-=iincx)
      {
        stemp += fabs(sx[i]);
      }
    }
  }

  return stemp;
} /* dasum_ */
