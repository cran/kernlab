#include <stdlib.h>
#include <string.h>
#include <R_ext/BLAS.h>
extern double mymin(double, double);
extern double mymax(double, double);
extern void *xmalloc(size_t);
/* LEVEL 1 BLAS */
/*extern double ddot_(int *, double *, int *, double *, int *);*/
/*extern int daxpy_(int *, double *, double *, int *, double *, int *);*/
/* LEVEL 2 BLAS */
/*extern int dsymv_(char *, int *, double *, double *, int *, double *, int *, double *, double *, int *);*/
/* MINPACK 2 */
extern void dbreakpt(int, double *, double *, double *, double *, int *, double *, double *);
extern void dgpstep(int, double *, double *, double *, double, double *, double *);

void dprsrch(int n, double *x, double *xl, double *xu, double *A, double *g, double *w)
{
/*
c     **********
c
c     Subroutine dprsrch
c
c     This subroutine uses a projected search to compute a step
c     that satisfies a sufficient decrease condition for the quadratic
c
c           q(s) = 0.5*s'*A*s + g'*s,
c
c     where A is a symmetric matrix and g is a vector. Given the 
c     parameter alpha, the step is
c
c           s[alpha] = P[x + alpha*w] - x,
c
c     where w is the search direction and P the projection onto the 
c     n-dimensional interval [xl,xu]. The final step s = s[alpha] 
c     satisfies the sufficient decrease condition
c
c           q(s) <= mu_0*(g'*s),
c
c     where mu_0 is a constant in (0,1).
c
c     The search direction w must be a descent direction for the
c     quadratic q at x such that the quadratic is decreasing
c     in the ray  x + alpha*w for 0 <= alpha <= 1.
c
c	parameters:
c
c       n is an integer variable.
c         On entry n is the number of variables.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x specifies the vector x.
c         On exit x is set to the final point P[x + alpha*w].
c
c       xl is a double precision array of dimension n.
c         On entry xl is the vector of lower bounds.
c         On exit xl is unchanged.
c
c       xu is a double precision array of dimension n.
c         On entry xu is the vector of upper bounds.
c         On exit xu is unchanged.
c
c       A is a double precision array of dimension n*n.
c         On entry A specifies the matrix A
c         On exit A is unchanged.
c
c       g is a double precision array of dimension n.
c         On entry g specifies the vector g.
c         On exit g is unchanged.
c
c       w is a double prevision array of dimension n.
c         On entry w specifies the search direction.
c         On exit w is the step s[alpha].
c
c     **********
*/

	double one = 1, zero = 0;

	/* Constant that defines sufficient decrease. */
	/* Interpolation factor. */
	double mu0 = 0.01, interpf = 0.5;
	
	double *wa1 = (double *) xmalloc(sizeof(double)*n);
	double *wa2 = (double *) xmalloc(sizeof(double)*n);

	/* Set the initial alpha = 1 because the quadratic function is 
	decreasing in the ray x + alpha*w for 0 <= alpha <= 1 */
	double alpha = 1, brptmin, brptmax, gts, q;
	int search = 1, nbrpt, nsteps = 0, i, inc = 1;	

	/* Find the smallest break-point on the ray x + alpha*w. */
	dbreakpt(n, x, xl, xu, w, &nbrpt, &brptmin, &brptmax);

	/* Reduce alpha until the sufficient decrease condition is
	satisfied or x + alpha*w is feasible. */
	while (search && alpha > brptmin)
	{

		/* Calculate P[x + alpha*w] - x and check the sufficient
		decrease condition. */
		nsteps++;
		dgpstep(n, x, xl, xu, alpha, w, wa1);
		F77_CALL(dsymv)("U", &n, &one, A, &n, wa1, &inc, &zero, wa2, &inc);
		gts = F77_CALL(ddot)(&n, g, &inc, wa1, &inc);
		q = 0.5*F77_CALL(ddot)(&n, wa1, &inc, wa2, &inc) + gts;
		if (q <= mu0*gts)
			search = 0;
		else
		
			/* This is a crude interpolation procedure that
			will be replaced in future versions of the code. */
			alpha *= interpf;
	}

	/* Force at least one more constraint to be added to the active
	set if alpha < brptmin and the full step is not successful. 
	There is sufficient decrease because the quadratic function 
	is decreasing in the ray x + alpha*w for 0 <= alpha <= 1. */
	if (alpha < 1 && alpha < brptmin)
		alpha = brptmin;

	/* Compute the final iterate and step. */
	dgpstep(n, x, xl, xu, alpha, w, wa1);
	F77_CALL(daxpy)(&n, &alpha, w, &inc, x, &inc);
	for (i=0;i<n;i++)
		x[i] = mymax(xl[i], mymin(x[i], xu[i]));
	memcpy(w, wa1, sizeof(double)*n);

	free(wa1);
	free(wa2);
}
