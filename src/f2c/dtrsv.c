#include "blas.h"

int dtrsv_(char *uplo, char *trans, char *diag, long *n, double *a,
           long *lda, double *x, long *incx)
{
  long info, i, j, ix, jx, kx;
  double temp;
  blasbool upper, notrans, nounit;

  /* pointers for testing */
  double *pa;

  /* Dereferenced input variables */
  long nn, dima, iincx;

  /* Dependencies */
  extern int xerbla_(char *, long *);

/*  Purpose   
    =======   

    DTRSV  solves one of the systems of equations   

       A*x = b,   or   A'*x = b,   

    where b and x are n element vectors and A is an n by n unit, or   
    non-unit, upper or lower triangular matrix.   

    No test for singularity or near-singularity is included in this   
    routine. Such tests must be performed before calling this routine.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the equations to be solved as   
             follows:   

                TRANS = 'N' or 'n'   A*x = b.   

                TRANS = 'T' or 't'   A'*x = b.   

                TRANS = 'C' or 'c'   A'*x = b.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
             triangular matrix and the strictly lower triangular part of 
             A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
             triangular matrix and the strictly upper triangular part of 
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of 
             A are not referenced either, but are assumed to be unity.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element right-hand side vector b. On exit, X is overwritten 
             with the solution vector x.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
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

  switch( *trans )
  {
    case 'N':
    case 'n':
      notrans = TRUE;
      break;
    case 'T':
    case 't':
    case 'C':
    case 'c':
      notrans = FALSE;
      break;
    default:
      notrans = TRUE;
      if( info == 0 )
        info = 2;
  }

  switch( *diag )
  {
    case 'N':
    case 'n':
      nounit = TRUE;
      break;
    case 'U':
    case 'u':
      nounit = FALSE;
      break;
    default:
      nounit = TRUE;
      if( info == 0 )
        info = 3;
  }

  if( info == 0 )
  {
    if (nn < 0) {
      info = 4;
    } else if (dima < MAX(1,nn)) {
      info = 6;
    } else if (iincx == 0) {
      info = 8;
    }
  }

  if (info != 0)
  {
    xerbla_("DTRSV ", &info);
    return 0;
  }

  if (nn != 0) /* Quick return if possible. */
  {
    /* Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (notrans) /* Form  x := inv( A )*x. */
    {
      if (upper)
      {
        if (iincx == 1)
        {
          for (pa=a+dima*(nn-1), j=nn-1; j>=0; j--, pa-=dima)
            if (x[j] != 0.0)
            {
              if (nounit)
                x[j] /= pa[j];
              temp = x[j];
              for (i = j - 1; i >= 0; i--)
                x[i] -= temp * pa[i];
            }
        }
        else
        {
          if (iincx >= 0) /* Set up the start point in X */
            jx = (nn - 1) * iincx;
          else
            jx = 0;
          for (pa=a+dima*(nn-1), j=nn-1; j>=0; j--, pa-=dima, jx-=iincx)
          {
            if (x[jx] != 0.0)
            {
              if (nounit)
                x[jx] /= pa[j];
              temp = x[jx];
              for (ix=jx, i=j-1; i>=0; i--)
              {
                ix -= iincx;
                x[ix] -= temp * pa[i];
              }
            }
          }
        }
      }
      else
      {
        if (iincx == 1)
        {
          for (pa=a, j=0; j<nn; j++, pa+=dima)
            if (x[j] != 0.0)
            {
              if (nounit)
                x[j] /= pa[j];
              temp = x[j];
              for (i = j + 1; i < nn; i++)
                x[i] -= temp * pa[i];
            }
        }
        else
        {
          if (iincx >= 0) /* Set up the start point in X */
            jx = 0;
          else
            jx = (1 - (nn)) * iincx;
          for (pa=a, j=0; j<nn; j++, pa+=dima, jx+=iincx)
          {
            if (x[jx] != 0.0)
            {
              if (nounit)
                x[jx] /= pa[j];
              temp = x[jx];
              for (ix=jx, i=j+1; i < nn; i++)
              {
                ix += iincx;
                x[ix] -= temp * pa[i];
              }
            }
          }
        }
      }
    }
    else /* Form  x := inv( A' )*x. */
    {
      if (upper)
      {
        if (iincx == 1)
        {
          for (pa=a, j=0; j<nn; j++, pa+=dima)
          {
            temp = x[j];
            for (i = 0; i < j; i++)
              temp -= pa[i] * x[i];
            if (nounit)
              temp /= pa[i];
            x[j] = temp;
          }
        }
        else
        {
          if (iincx >= 0) /* Set up the start point in X */
            kx = 0;
          else
            kx = (1 - nn) * iincx;
          for (pa=a, jx=kx, j=0; j<nn; j++, pa+=dima, jx+=iincx)
          {
            temp = x[jx];
            for (ix=kx, i=0; i<j; i++, ix+=iincx)
              temp -= pa[i] * x[ix];
            if (nounit)
              temp /= pa[i];
            x[jx] = temp;
          }
        }
      }
      else
      {
        if (iincx == 1)
        {
          for (pa=a+dima*(nn-1), j=nn-1; j>=0; j--, pa-=dima)
          {
            temp = x[j];
            for (i = nn-1; i > j; i--)
              temp -= pa[i] * x[i];
            if (nounit)
              temp /= pa[i];
            x[j] = temp;
          }
        }
        else
        {
          if (iincx >= 0) /* Set up the start point in X */
            kx = (nn - 1) * iincx;
          else
            kx = 0;
          for (pa=a+dima*(nn-1), jx=kx, j=nn-1; j>=0; j--, pa-=dima, jx-=iincx)
          {
            temp = x[jx];
            for (ix=kx, i=nn-1; i>j; i--, ix-=iincx)
              temp -= pa[i] * x[ix];
            if (nounit)
              temp /= pa[i];
            x[jx] = temp;
          }
        }
      }
    }
  }

  return 0;
} /* dtrsv_ */
