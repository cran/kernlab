#include <stdlib.h>
#include <math.h>
#include <string.h>

extern void *xmalloc(size_t);
/* LEVEL 1 BLAS */
extern int daxpy_(int *, double *, double *, int *, double *, int *);
extern double ddot_(int *, double *, int *, double *, int *);
extern double dnrm2_(int *, double *, int *);
extern int dscal_(int *, double *, double *, int *);
/* LEVEL 2 BLAS */
extern int dtrsv_(char *, char *, char *, int *, double *, int *, double *,  int *);
extern int dsymv_(char *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
/* MINPACK 2 */
extern void dtrqsol(int, double *, double *, double , double *);

void dtrpcg(int n, double *A, double *g, double delta, double *L, double tol, double stol, double *w, int *iters, int *info)
{
/*
c     *********
c
c     Subroutine dtrpcg
c
c     Given a dense symmetric positive semidefinite matrix A, this
c     subroutine uses a preconditioned conjugate gradient method to find
c     an approximate minimizer of the trust region subproblem
c
c           min { q(s) : || L'*s || <= delta }.
c
c     where q is the quadratic
c
c           q(s) = 0.5*s'*A*s + g'*s,
c
c     This subroutine generates the conjugate gradient iterates for
c     the equivalent problem
c
c           min { Q(w) : || w || <= delta }.
c
c     where Q is the quadratic defined by
c
c           Q(w) = q(s),      w = L'*s.
c
c     Termination occurs if the conjugate gradient iterates leave
c     the trust region, a negative curvature direction is generated,
c     or one of the following two convergence tests is satisfied.
c
c     Convergence in the original variables:
c
c           || grad q(s) || <= tol
c
c     Convergence in the scaled variables:
c
c           || grad Q(w) || <= stol
c
c     Note that if w = L'*s, then L*grad Q(w) = grad q(s).
c
c     parameters:
c
c       n is an integer variable.
c         On entry n is the number of variables.
c         On exit n is unchanged.
c
c       A is a double precision array of dimension n*n.
c         On entry A specifies the matrix A.
c         On exit A is unchanged.
c
c       g is a double precision array of dimension n.
c         On entry g must contain the vector g.
c         On exit g is unchanged.
c
c       delta is a double precision variable.
c         On entry delta is the trust region size.
c         On exit delta is unchanged.
c
c       L is a double precision array of dimension n*n.
c         On entry L need not to be specified.
c         On exit the lower triangular part of L contains the matrix L.
c
c       tol is a double precision variable.
c         On entry tol specifies the convergence test
c            in the un-scaled variables.
c         On exit tol is unchanged
c
c       stol is a double precision variable.
c         On entry stol specifies the convergence test
c            in the scaled variables.
c         On exit stol is unchanged 
c
c       w is a double precision array of dimension n.
c         On entry w need not be specified.
c         On exit w contains the final conjugate gradient iterate.
c
c       iters is an integer variable.
c         On entry iters need not be specified.
c         On exit iters is set to the number of conjugate
c            gradient iterations.
c
c       info is an integer variable.
c         On entry info need not be specified.
c         On exit info is set as follows:
c
c             info = 1  Convergence in the original variables.
c                       || grad q(s) || <= tol
c
c             info = 2  Convergence in the scaled variables.
c                       || grad Q(w) || <= stol
c
c             info = 3  Negative curvature direction generated.
c                       In this case || w || = delta and a direction
c
c                       of negative curvature w can be recovered by
c                       solving L'*w = p.
c
c             info = 4  Conjugate gradient iterates exit the
c                       trust region. In this case || w || = delta.
c
c             info = 5  Failure to converge within itermax(n) iterations.
c
c     **********
*/
	int i, inc = 1;
	double one = 1, zero = 0, alpha, malpha, beta, ptq, rho;
	double *p, *q, *t, *r, *z, sigma, rtr, rnorm, rnorm0, tnorm;
	p = (double *) xmalloc(sizeof(double)*n);
	q = (double *) xmalloc(sizeof(double)*n);
	t = (double *) xmalloc(sizeof(double)*n);
	r = (double *) xmalloc(sizeof(double)*n);
	z = (double *) xmalloc(sizeof(double)*n);

	/* Initialize the iterate w and the residual r.
	Initialize the residual t of grad q to -g.
	Initialize the residual r of grad Q by solving L*r = -g.
	Note that t = L*r. */
	for (i=0;i<n;i++)
	{
		w[i] = 0;
		r[i] = t[i] = -g[i];
	}
	dtrsv_("L", "N", "N", &n, L, &n, r, &inc);

	/* Initialize the direction p. */
	memcpy(p, r, sizeof(double)*n);
	
	/* Initialize rho and the norms of r and t. */
	rho = ddot_(&n, r, &inc, r, &inc);
	rnorm0 = sqrt(rho);

	/* Exit if g = 0. */
	*iters = 0;
	if (rnorm0 <= stol)
	{
		*info = 1;
		goto return0;
	}

	for (*iters=1;*iters<=n;*iters++)
	{
		
		/* Compute z by solving L'*z = p. */
		memcpy(z, p, sizeof(double)*n);
		dtrsv_("L", "T", "N", &n, L, &n, z, &inc);

		/* Compute q by solving L*q = A*z and save L*q for
		use in updating the residual t.	*/
		dsymv_("U", &n, &one, A, &n, z, &inc, &zero, q, &inc);
		memcpy(z, q, sizeof(double)*n);
		dtrsv_("L", "N", "N", &n, L, &n, q, &inc);
		
		/* Compute alpha and determine sigma such that the trust region
		constraint || w + sigma*p || = delta is satisfied. */
		ptq = ddot_(&n, p, &inc, q, &inc);
		if (ptq > 0)
			alpha = rho/ptq;
		else
			alpha = 0;
		dtrqsol(n, w, p, delta, &sigma);
	
		/* Exit if there is negative curvature or if the
		iterates exit the trust region.	*/
		if (ptq <= 0 || alpha >= sigma)
		{
			daxpy_(&n, &sigma, p, &inc, w, &inc);
			if (ptq <= 0)
				*info = 3;
			else
				*info = 4;
			goto return0;
		}

		/* Update w and the residuals r and t.
		Note that t = L*r. */
		malpha = -alpha;
		daxpy_(&n, &alpha, p, &inc, w, &inc);
		daxpy_(&n, &malpha, q, &inc, r, &inc);
		daxpy_(&n, &malpha, z, &inc, t,&inc);
	
		/* Exit if the residual convergence test is satisfied. */
		rtr = ddot_(&n, r, &inc, r, &inc);
		rnorm = sqrt(rtr);
		tnorm = sqrt(ddot_(&n, t, &inc, t, &inc));
		if (tnorm <= tol)
		{
			*info = 1;
			goto return0;
		}
		if (rnorm <= stol)
		{
			*info = 2;
			goto return0;
		}

		/* Compute p = r + beta*p and update rho. */
		beta = rtr/rho;
		dscal_(&n, &beta, p, &inc);
		daxpy_(&n, &one, r, &inc, p, &inc);
		rho = rtr;
	}

	/* iters > itermax = n */
	*info = 5;
return0:
	free(p);
	free(q);
	free(r);
	free(t);
	free(z);
} 
