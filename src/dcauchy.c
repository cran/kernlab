#include <stdlib.h>
#include <R_ext/BLAS.h>

extern void *xmalloc(size_t);
/* LEVEL 1 BLAS */
/* extern double ddot_(int *, double *, int *, double *, int *);
 extern double dnrm2_(int *, double *, int *); */
/* LEVEL 2 BLAS */
/* extern int dsymv_(char *, int *, double *, double *, int *, double *, int *, double *, double *, int *); */
/* MINPACK 2 */
extern void dbreakpt(int, double *, double *, double *, double *, int *, double *, double *);
extern void dgpstep(int, double *, double *, double *, double, double *, double *);

void dcauchy(int n, double *x, double *xl, double *xu, double *A, double *g, double delta, double *alpha, double *s)
{
/*
c     **********
c
c     Subroutine dcauchy
c
c     This subroutine computes a Cauchy step that satisfies a trust
c     region constraint and a sufficient decrease condition.
c
c     The Cauchy step is computed for the quadratic
c
c           q(s) = 0.5*s'*A*s + g'*s,
c
c     where A is a symmetric matrix , and g is a vector. Given a 
c     parameter alpha, the Cauchy step is
c
c           s[alpha] = P[x - alpha*g] - x,
c
c     with P the projection onto the n-dimensional interval [xl,xu].
c     The Cauchy step satisfies the trust region constraint and the
c     sufficient decrease condition
c
c           || s || <= delta,      q(s) <= mu_0*(g'*s),
c
c     where mu_0 is a constant in (0,1).
c
c	parameters:
c
c       n is an integer variable.
c         On entry n is the number of variables.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x specifies the vector x.
c         On exit x is unchanged.
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
c         On entry A specifies the matrix A.
c         On exit A is unchanged.
c
c       g is a double precision array of dimension n.
c         On entry g specifies the gradient g.
c         On exit g is unchanged.
c
c       delta is a double precision variable.
c         On entry delta is the trust region size.
c         On exit delta is unchanged.
c
c       alpha is a double precision variable.
c         On entry alpha is the current estimate of the step.
c         On exit alpha defines the Cauchy step s[alpha].
c
c       s is a double precision array of dimension n.
c         On entry s need not be specified.
c         On exit s is the Cauchy step s[alpha].
c
c     **********
*/

	double one = 1, zero = 0;

	/* Constant that defines sufficient decrease.
	Interpolation and extrapolation factors. */
	double mu0 = 0.01, interpf = 0.1, extrapf = 10;
	
	int search, interp, nbrpt, nsteps = 1, i, inc = 1;
	double alphas, brptmax, brptmin, gts, q; 
	double *wa = (double *) xmalloc(sizeof(double)*n);
	
	/* Find the minimal and maximal break-point on x - alpha*g. */
	for (i=0;i<n;i++)
		wa[i] = -g[i];
	dbreakpt(n, x, xl, xu, wa, &nbrpt, &brptmin, &brptmax);
	
	/* Evaluate the initial alpha and decide if the algorithm
	must interpolate or extrapolate. */
	dgpstep(n, x, xl, xu, -(*alpha), g, s);
	if (F77_CALL(dnrm2)(&n, s, &inc) > delta)
		interp = 1;
	else
	{
		F77_CALL(dsymv)("U", &n, &one, A, &n, s, &inc, &zero, wa, &inc);
		gts = F77_CALL(ddot)(&n, g, &inc, s, &inc);
		q = 0.5*F77_CALL(ddot)(&n, s, &inc, wa, &inc) + gts;
		interp = q >= mu0*gts ? 1 : 0;
	}
	
	/* Either interpolate or extrapolate to find a successful step. */
	if (interp)
	{

		/* Reduce alpha until a successful step is found. */
		search = 1;
		while (search)
		{

			/* This is a crude interpolation procedure that
			will be replaced in future versions of the code. */
			nsteps++;
			(*alpha) *= interpf;
			dgpstep(n, x, xl, xu, -(*alpha), g, s);
			if (F77_CALL(dnrm2)(&n, s, &inc) <= delta)
			{
				F77_CALL(dsymv)("U", &n, &one, A, &n, s, &inc, &zero, wa, &inc);
				gts = F77_CALL(ddot)(&n, g, &inc, s, &inc);
				q = 0.5 * F77_CALL(ddot)(&n, s, &inc, wa, &inc) + gts;
				search = q > mu0*gts ? 1 : 0;
			} 			
		}	
	}
	else
	{
		search = 1;
		alphas = *alpha;
	
		/* Increase alpha until a successful step is found. */
		while (search && (*alpha) <= brptmax)
		{
		
			/* This is a crude extrapolation procedure that
			will be replaced in future versions of the code. */
			nsteps++;
			alphas = *alpha;
			(*alpha) *= extrapf;
			dgpstep(n, x, xl, xu, -(*alpha), g, s);
			if (F77_CALL(dnrm2)(&n, s, &inc) <= delta)
			{
				F77_CALL(dsymv)("U", &n, &one, A, &n, s, &inc, &zero, wa, &inc);
				gts = F77_CALL(ddot)(&n, g, &inc, s, &inc);
				q = 0.5 * F77_CALL(ddot)(&n, s, &inc, wa, &inc) + gts;
				search = q < mu0*gts ? 1 : 0;
			}
			else
				search = 0;					
		}
		*alpha = alphas;
		dgpstep(n, x, xl, xu, -(*alpha), g, s);
	}

	free(wa);
}
