#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern void *xmalloc(size_t);
extern double mymin(double, double);
extern double mymax(double, double);
extern int ufv(int, double *, double *);
extern int ugrad(int, double *, double *);
extern int uhes(int, double *, double **);
/* LEVEL 1 BLAS */
extern double dnrm2_(int *, double *, int *);
extern double ddot_(int *, double *, int *, double *, int *);
/* LEVEL 2 BLAS */
extern int dsymv_(char *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
/* MINPACK 2 */
extern double dgpnrm(int, double *, double *, double *, double *);
extern void dcauchy(int, double *, double *, double *, double *, double *, double, double *, double *, double *);
extern void dspcg(int, double *, double *, double *, double *, double *, double, double, double *, int *);

void dtron(int n, double *x, double *xl, double *xu, double gtol, double frtol, double fatol, double fmin, int maxfev, double cgtol) 
{
/*
c     *********
c
c     Subroutine dtron
c
c     The optimization problem of BSVM is a bound-constrained quadratic
c     optimization problem and its Hessian matrix is positive semidefinite. 
c     We modified the optimization solver TRON by Chih-Jen Lin and
c     Jorge More' into this version which is suitable for this
c     special case.
c
c     This subroutine implements a trust region Newton method for the
c     solution of large bound-constrained quadratic optimization problems
c
c           min { f(x)=0.5*x'*A*x + g0'*x : xl <= x <= xu }
c
c     where the Hessian matrix A is dense and positive semidefinite. The
c     user must define functions which evaluate the function, gradient, 
c     and the Hessian matrix.
c
c     The user must choose an initial approximation x to the minimizer,
c     lower bounds, upper bounds, quadratic terms, linear terms, and
c     constants about termination criterion.
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
c       gtol is a double precision variable.
c         On entry gtol specifies the relative error of the projected 
c            gradient.
c         On exit gtol is unchanged. 
c
c       frtol is a double precision variable.
c         On entry frtol specifies the relative error desired in the
c            function. Convergence occurs if the estimate of the
c            relative error between f(x) and f(xsol), where xsol
c            is a local minimizer, is less than frtol.
c         On exit frtol is unchanged.
c
c       fatol is a double precision variable.
c         On entry fatol specifies the absolute error desired in the
c            function. Convergence occurs if the estimate of the
c            absolute error between f(x) and f(xsol), where xsol
c            is a local minimizer, is less than fatol.
c         On exit fatol is unchanged.
c
c       fmin is a double precision variable.
c         On entry fmin specifies a lower bound for the function.
c            The subroutine exits with a warning if f < fmin.
c         On exit fmin is unchanged.
c
c       maxfev is an integer variable.
c         On entry maxfev specifies the limit of function evaluations.
c         On exit maxfev is unchanged.
c
c       cgtol is a double precision variable.
c         On entry gqttol specifies the convergence criteria for
c            subproblems.
c         On exit gqttol is unchanged.
c
c     **********
*/

	/* Parameters for updating the iterates. */
	double eta0 = 1e-4, eta1 = 0.25, eta2 = 0.75;
	
	/* Parameters for updating the trust region size delta. */
	double sigma1 = 0.25, sigma2 = 0.5, sigma3 = 4;

	double p5 = 0.5, one = 1;
	double gnorm, gnorm0, delta, snorm;
	double alphac = 1, alpha, f, fc, prered, actred, gs;
	int search = 1, iter = 1, info, inc = 1;	
	double *xc = (double *) xmalloc(sizeof(double)*n);
	double *s = (double *) xmalloc(sizeof(double)*n);
	double *wa = (double *) xmalloc(sizeof(double)*n); 	
	double *g = (double *) xmalloc(sizeof(double)*n);
	double *A = NULL;

	uhes(n, x, &A);
	ugrad(n, x, g);	
	ufv(n, x, &f);	
	gnorm0 = dnrm2_(&n, g, &inc);
	delta = 1000*gnorm0;
	gnorm = dgpnrm(n, x, xl, xu, g);
	if (gnorm <= gtol*gnorm0)
	{
		/*
		printf("CONVERGENCE: GTOL TEST SATISFIED\n");
		*/
		search = 0;
	}

	while (search)
	{

		/* Save the best function value and the best x. */   
		fc = f;
		memcpy(xc, x, sizeof(double)*n);
		
		/* Compute the Cauchy step and store in s. */		
		dcauchy(n, x, xl, xu, A, g, delta, &alphac, s, wa);
		
		/* Compute the projected Newton step. */		
		dspcg(n, x, xl, xu, A, g, delta, cgtol, s, &info);
		if (ufv(n, x, &f) > maxfev)
		{
			/*
			printf("ERROR: NFEV > MAXFEV\n");
			*/
			search = 0;
			continue;
		}

		/* Compute the predicted reduction. */
		memcpy(wa, g, sizeof(double)*n);
		dsymv_("U", &n, &p5, A, &n, s, &inc, &one, wa, &inc);
		prered = -ddot_(&n, s, &inc, wa, &inc);
                        
		/* Compute the actual reduction. */
	        actred = fc - f;                
                                    
		/* On the first iteration, adjust the initial step bound. */
		snorm = dnrm2_(&n, s, &inc);
		if (iter == 1)
			delta = mymin(delta, snorm);

		/* Compute prediction alpha*snorm of the step. */
		gs = ddot_(&n, g, &inc, s, &inc);

		if (f - fc - gs <= 0)
			alpha = sigma3;
		else
			alpha = mymax(sigma1, -0.5*(gs/(f - fc - gs)));

		/* Update the trust region bound according to the ratio
		of actual to predicted reduction. */
		if (actred < eta0*prered)
		
			/* Reduce delta. Step is not successful. */
			delta = mymin(mymax(alpha, sigma1)*snorm, sigma2*delta);
		else 
		{
			if (actred < eta1*prered)

				/* Reduce delta. Step is not sufficiently successful. */
				delta = mymax(sigma1*delta, mymin(alpha*snorm, sigma2*delta));
			else 
				if (actred < eta2*prered)
				
					/* The ratio of actual to predicted reduction is in
					the interval (eta1,eta2). We are allowed to either
					increase or decrease delta. */
					delta = mymax(sigma1*delta, mymin(alpha*snorm, sigma3*delta));
				else
				
					/* The ratio of actual to predicted reduction exceeds eta2.
					Do not decrease delta. */
					delta = mymax(delta, mymin(alpha*snorm, sigma3*delta));
		}

		/* Update the iterate. */
		if (actred > eta0*prered) 
		{
		
			/* Successful iterate. */
			iter++;
			/*
			uhes(n, x, &A);
			*/
			ugrad(n, x, g);
			gnorm = dgpnrm(n, x, xl, xu, g);		
			if (gnorm <= gtol*gnorm0)
        		{
				/*
				printf("CONVERGENCE: GTOL = %g TEST SATISFIED\n", gnorm/gnorm0);
				*/
				search = 0;
		                continue;
			}
		}
		else
		{

			/* Unsuccessful iterate. */
			memcpy(x, xc, sizeof(double)*n);
			f = fc;
		}
	
		/* Test for convergence */
		if (f < fmin)
		{
			printf("WARNING: F .LT. FMIN\n");
			search = 0; /* warning */
			continue;
		}
		if (fabs(actred) <= fatol && prered <= fatol)
		{
			printf("CONVERGENCE: FATOL TEST SATISFIED\n");
			search = 0;
			continue;
		}
		if (fabs(actred) <= frtol*fabs(f) && prered <= frtol*fabs(f))
		{
			/*
			printf("CONVERGENCE: FRTOL TEST SATISFIED\n");		
			*/
			search = 0;
			continue;
		}
	}	

	free(g);
	free(xc);
	free(s);
	free(wa);
}
