#include "blas.h"

int dsymv_(char *uplo, long *n, double *alpha, double *a, long *lda,
           double *x, long *incx, double *beta, double *y, long *incy)
{
  long info, i, j, ix, iy, jx, jy, kx, ky;
  double temp1, temp2;
  blasbool upper;

  /* pointers for testing */
  double *pa;

  /* Dereferenced input variables */
  long nn, dima, iincx, iincy;
  double aalpha, bbeta;

  /* Dependencies */
  extern int xerbla_(char *, long *);

/*  Purpose   
    =======   

    DSYMV  performs the matrix-vector  operation   

       y := alpha*A*x + beta*y,   

    where alpha and beta are scalars, x and y are n element vectors and   
    A is an n by n symmetric matrix.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
             triangular part of the symmetric matrix and the strictly   
             lower triangular part of A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
             triangular part of the symmetric matrix and the strictly   
             upper triangular part of A is not referenced.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y. On exit, Y is overwritten by the updated   
             vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   
*/

  /* Dereference the inputs */
  nn = *n;
  dima = *lda;
  iincx = *incx;
  iincy = *incy;
  aalpha = *alpha;
  bbeta = *beta;

  info = 0;

  switch( *uplo )
  {
    case 'L':
    case 'l':
      upper = FALSE;
      break;
    case 'U':
    case 'u':
      upper = TRUE;
      break;
    default:
      upper = FALSE;
      info = 1;
  }

  if (info == 0)
  {
    if (nn < 0) {
      info = 2;
    } else if (dima < MAX(1,nn)) {
      info = 5;
    } else if (iincx == 0) {
      info = 7;
    } else if (iincy == 0) {
      info = 10;
    }
  }

  if (info != 0)
  {
    xerbla_("DSYMV ", &info);
    return 0;
  }

  /* Quick return if possible. */

  if (nn != 0 && (aalpha != 0.0 || bbeta != 1.0))
  {

    /* Set up the start points in  X  and  Y. */

    if (iincx > 0)
      kx = 0;
    else
      kx = (1 - nn) * iincx;
    if (iincy > 0)
      ky = 0;
    else
      ky = (1 - nn) * iincy;

    /* Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A. */

    /* First form  y := beta*y. */

    if (bbeta != 1.0)
    {
      if (iincy == 1)
      {
        if (bbeta == 0.0)
        {
          for (i = 0; i < nn; i++)
            y[i] = 0.0;
        }
        else
        {
          for (i = 0; i < nn; i++)
            y[i] = bbeta * y[i];
        }
      }
      else
      {
        iy = ky;
        if (bbeta == 0.0)
        {
          for (i=0; i<nn; i++, iy+=iincy)
            y[iy] = 0.0;
        }
        else
        {
          for (i=0; i<nn; i++, iy+=iincy)
            y[iy] = bbeta * y[iy];
        }
      }
    }

    if (aalpha != 0.0)
    {
      if (upper) /* Form  y  when A is stored in upper triangle. */
      {
        if (iincx == 1 && iincy == 1)
        {
          for (pa=a, j=0; j<nn; j++, pa+=dima)
          {
            temp1 = aalpha * x[j];
            temp2 = 0.0;
            for (i = 0; i < j; i++)
            {
              y[i] += temp1 * pa[i];
              temp2 += pa[i] * x[i];
            }
            y[i] += temp1 * pa[i] + aalpha * temp2;
          }
        }
        else
        {
          for (pa=a, jx=kx, jy=ky, j=0; j<nn; j++, pa+=dima, jx+=iincx, jy+=iincy)
          {
            temp1 = aalpha * x[jx];
            temp2 = 0.0;
            for (ix=kx, iy=ky, i=0; i<j; i++, ix+=iincx, iy+=iincy)
            {
              y[iy] += temp1 * pa[i];
              temp2 += pa[i] * x[ix];
            }
            y[jy] += temp1 * pa[j] + aalpha * temp2;   /* ??? diff indices? */
          }
        }
      }
      else /* Form  y  when A is stored in lower triangle. */
      {
        if (iincx == 1 && iincy == 1)
        {
          for (pa=a, j=0; j<nn; j++, pa+=dima)
          {
            temp1 = aalpha * x[j];
            temp2 = 0.0;
            y[j] += temp1 * pa[j];
            for (i=j+1; i<nn; i++)
            {
              y[i] += temp1 * pa[i];
              temp2 += pa[i] * x[i];
            }
            y[j] += aalpha * temp2;
          }
        }
        else
        {
          for (pa=a, jx=kx, jy=ky, j=0; j<nn; j++, jx+=iincx, jy+=iincy, pa+=dima)
          {
            temp1 = aalpha * x[jx];
            temp2 = 0.0;
            y[jy] += temp1 * pa[j];
            for (ix=jx, iy=jy, i=j+1; i<nn; i++)
            {
              ix += iincx;
              iy += iincy;
              y[iy] += temp1 * pa[i];
              temp2 += pa[i] * x[ix];
            }
            y[jy] += aalpha * temp2;
          }
        }
      }
    }
  }

  return 0;
} /* dsymv_ */
