#include <math.h>
#include <R_ext/BLAS.h>

extern double mymax(double, double);
/* LEVEL 1 BLAS */
/*extern double ddot_(int *, double *, int *, double *, int *);*/

void dtrqsol(int n, double *x, double *p, double delta, double *sigma)
{
/*
c     **********
c
c     Subroutine dtrqsol
c
c     This subroutine computes the largest (non-negative) solution
c     of the quadratic trust region equation
c
c           ||x + sigma*p|| = delta.
c
c     The code is only guaranteed to produce a non-negative solution
c     if ||x|| <= delta, and p != 0. If the trust region equation has
c     no solution, sigma = 0.
c
c	parameters:
c
c       n is an integer variable.
c         On entry n is the number of variables.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x must contain the vector x.
c         On exit x is unchanged.
c
c       p is a double precision array of dimension n.
c         On entry p must contain the vector p.
c         On exit p is unchanged.
c
c       delta is a double precision variable.
c         On entry delta specifies the scalar delta.
c         On exit delta is unchanged.
c
c       sigma is a double precision variable.
c         On entry sigma need not be specified.
c         On exit sigma contains the non-negative solution.
c
c     **********
*/
	int inc = 1;
	double dsq = delta*delta, ptp, ptx, rad, xtx;
	ptx = F77_CALL(ddot)(&n, p, &inc, x, &inc);
	ptp = F77_CALL(ddot)(&n, p, &inc, p, &inc);
	xtx = F77_CALL(ddot)(&n, x, &inc, x, &inc);

	/* Guard against abnormal cases. */
	rad = ptx*ptx + ptp*(dsq - xtx);
	rad = sqrt(mymax(rad, 0));
	if (ptx > 0)
		*sigma = (dsq - xtx)/(ptx + rad);
	else
		if (rad > 0)
			*sigma = (rad - ptx)/ptp;
		else
			*sigma = 0;
}
