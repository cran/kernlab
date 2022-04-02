#include <stdlib.h>
#ifndef USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <R_ext/BLAS.h>
extern void *xmalloc(size_t);
extern double mymin(double, double);
extern double mymax(double, double);
/* LEVEL 1 BLAS */
/*extern double dnrm2_(int *, double *, int *);*/
/* LEVEL 2 BLAS */
/*extern int dsymv_(char *, int *, double *, double *, int *, double *, int *, double *, double *, int *);*/
/*extern void dtrsv_(char *, char *, char *, int *, double *, int *, double *, int *);*/
/* MINPACK 2 */
extern void dprsrch(int, double *, double *, double *, double *, double *, double *);
extern double dprecond(int, double *, double *);
extern void dtrpcg(int, double*, double *, double, double *, double, double, double *, int *, int *);

void dspcg(int n, double *x, double *xl, double *xu, double *A, double *g, double delta, double rtol, double *s, int *info)
{
/*
c     *********
c
c     Subroutine dspcg
c
c     This subroutine generates a sequence of approximate minimizers
c     for the subproblem
c
c           min { q(x) : xl <= x <= xu }.
c
c     The quadratic is defined by
c
c           q(x[0]+s) = 0.5*s'*A*s + g'*s,
c
c     where x[0] is a base point provided by the user, A is a symmetric
c     positive semidefinite dense matrix, and g is a vector.
c
c     At each stage we have an approximate minimizer x[k], and generate
c     a direction p[k] by solving the subproblem
c
c           min { q(x[k]+p) : || p || <= delta, s(fixed) = 0 },
c
c     where fixed is the set of variables fixed at x[k], delta is the
c     trust region bound.
c
c           B = A(free:free),
c
c     where free is the set of free variables at x[k]. Given p[k],
c     the next minimizer x[k+1] is generated by a projected search.
c
c     The starting point for this subroutine is x[1] = x[0] + s, where
c     x[0] is a base point and s is the Cauchy step.
c
c     The subroutine converges when the step s satisfies
c
c           || (g + A*s)[free] || <= rtol*|| g[free] ||
c
c     In this case the final x is an approximate minimizer in the
c     face defined by the free variables.
c
c     The subroutine terminates when the trust region bound does
c     not allow further progress, that is, || L'*p[k] || = delta.
c     In this case the final x satisfies q(x) < q(x[k]).
c
c	parameters:
c
c       n is an integer variable.
c         On entry n is the number of variables.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x specifies the vector x.
c         On exit x is the final minimizer.
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
c         On entry g must contain the vector g.
c         On exit g is unchanged.
c
c       delta is a double precision variable.
c         On entry delta is the trust region size.
c         On exit delta is unchanged.
c
c       rtol is a double precision variable.
c         On entry rtol specifies the accuracy of the final minimizer.
c         On exit rtol is unchanged.
c
c       s is a double precision array of dimension n.
c         On entry s is the Cauchy step.
c         On exit s contain the final step.
c
c       info is an integer variable.
c         On entry info need not be specified.
c         On exit info is set as follows:
c
c             info = 1  Convergence. The final step s satisfies
c                       || (g + A*s)[free] || <= rtol*|| g[free] ||,
c                       and the final x is an approximate minimizer
c                       in the face defined by the free variables.
c
c             info = 2  Termination. The trust region bound does
c                       not allow further progress.
*/
	int i, j, nfaces, nfree, inc = 1, infotr, iters = 0, itertr;
	double gfnorm, gfnormf, stol = 1e-16, alpha;
	double one = 1, zero = 0;
	
	double *B = (double *) xmalloc(sizeof(double)*n*n);
	double *L = (double *) xmalloc(sizeof(double)*n*n);
	double *w = (double *) xmalloc(sizeof(double)*n);
	double *wa = (double *) xmalloc(sizeof(double)*n);
	double *wxl = (double *) xmalloc(sizeof(double)*n);
	double *wxu = (double *) xmalloc(sizeof(double)*n);	
	int *indfree = (int *) xmalloc(sizeof(int)*n);
	double *gfree = (double *) xmalloc(sizeof(double)*n);

	/* Compute A*(x[1] - x[0]) and store in w. */
	F77_CALL(dsymv)("U", &n, &one, A, &n, s, &inc, &zero, w, &inc FCONE);
      
	/* Compute the Cauchy point. */
	for (j=0;j<n;j++)
	{
		x[j] += s[j];
		x[j] = mymax(xl[j], mymin(x[j], xu[j]));
	}      
                 
	/* Start the main iteration loop.
	There are at most n iterations because at each iteration
	at least one variable becomes active. */
	for (nfaces=0;nfaces<n;nfaces++)
	{
	
		/* Determine the free variables at the current minimizer.
		The indices of the free variables are stored in the first
		n free positions of the array indfree. */
		nfree = 0;
		for (j=0;j<n;j++)
			if (xl[j] < x[j] && x[j] < xu[j])
				indfree[nfree++] = j;
			
		/* Exit if there are no free constraints. */
		if (nfree == 0)
		{
			*info = 1;
			goto return0;
		}

		/* Obtain the submatrix of A for the free variables. 
		Compute the gradient grad q(x[k]) = g + A*(x[k] - x[0]),
		of q at x[k] for the free variables.
		Recall that w contains  A*(x[k] - x[0]).
		Compute the norm of the reduced gradient Z'*g. */
		for (i=0;i<nfree;i++)
		{
			for (j=0;j<nfree;j++)
				B[i*nfree+j] = A[indfree[i]*n+indfree[j]];
			wa[i] = g[indfree[i]];
			gfree[i] = w[indfree[i]] + wa[i];
		}
		gfnorm = F77_CALL(dnrm2)(&nfree, wa, &inc);

		alpha = dprecond(nfree, B, L);
		dtrpcg(nfree, B, gfree, delta, L, rtol*gfnorm, stol, w, &itertr, &infotr);
		iters += itertr;
		F77_CALL(dtrsv)("L", "T", "N", &nfree, L, &nfree, w, &inc FCONE FCONE FCONE);

		/* Use a projected search to obtain the next iterate.
		The projected search algorithm stores s[k] in w. */
		for (j=0;j<nfree;j++)
		{
			wa[j] = x[indfree[j]];
			wxl[j] = xl[indfree[j]];
			wxu[j] = xu[indfree[j]];
		}
		dprsrch(nfree, wa, wxl, wxu, B, gfree, w);
		
		/* Update the minimizer and the step.
		Note that s now contains x[k+1] - x[0].	*/
		for (j=0;j<nfree;j++)
		{
			x[indfree[j]] = wa[j];
			s[indfree[j]] += w[j];
		}

		/* Compute A*(x[k+1] - x[0]) and store in w. */
		F77_CALL(dsymv)("U", &n, &one, A, &n, s, &inc, &zero, w, &inc FCONE);
         
		/* Compute the gradient grad q(x[k+1]) = g + A*(x[k+1] - x[0])
		of q at x[k+1] for the free variables. */
		for (j=0;j<nfree;j++)
			gfree[j] = w[indfree[j]] + g[indfree[j]];
		gfnormf = F77_CALL(dnrm2)(&nfree, gfree, &inc);

		/* Convergence and termination test.
		We terminate if the preconditioned conjugate gradient method
		encounters a direction of negative curvature, or
		if the step is at the trust region bound. */
		if (gfnormf <= rtol*gfnorm)
		{
			*info = 1;
			goto return0;
		}
		else
			if (infotr == 3 || infotr == 4)
			{
				*info = 2;
				goto return0;
			}
	}

return0:
	free(B);
	free(L);
	free(w);
	free(wa);
	free(wxl);
	free(wxu);
	free(indfree);
	free(gfree);		
} 
