#include <stdlib.h>
#include <string.h>
#include <R_ext/BLAS.h>
/* LEVEL 1 BLAS */
/*extern double ddot_(int *, double *, int *, double *, int *); */
/* LEVEL 2 BLAS */
/*extern int dsymv_(char *, int *, double *, double *, int *, double *, int *, double *, double *, int *);*/
/* MINPACK 2 */
extern void dtron(int, double *, double *, double *, double, double, double, double, int, double);

struct BQP
{
        double eps;
        int n;
        double *x, *C, *Q, *p;
};

int nfev, inc = 1;
double one = 1, zero = 0, *A, *g0;

int uhes(int n, double *x, double **H)
{
	*H = A;
	return 0;
}
int ugrad(int n, double *x, double *g)
{
	/* evaluate the gradient g = A*x + g0 */
	memcpy(g, g0, sizeof(double)*n);
	F77_CALL(dsymv)("U", &n, &one, A, &n, x, &inc, &one, g, &inc);
	return 0;
}
int ufv(int n, double *x, double *f)
{
	/* evaluate the function value f(x) = 0.5*x'*A*x + g0'*x */  
	double *t = (double *) malloc(sizeof(double)*n);
	F77_CALL(dsymv)("U", &n, &one, A, &n, x, &inc, &zero, t, &inc);
	*f = F77_CALL(ddot)(&n, x, &inc, g0, &inc) + 0.5 * F77_CALL(ddot)(&n, x, &inc, t, &inc);
	free(t);
	return ++nfev;
}

void solvebqp(struct BQP *qp)
{
	/* driver for positive semidefinite quadratic programing version
	of tron */
	int i, n, maxfev;
	double *x, *xl, *xu;
	double frtol, fatol, fmin, gtol, cgtol;

	n = qp->n;
	maxfev = 1000; /* ? */
	nfev = 0;

	x = qp->x;
	xu = qp->C;
	A = qp->Q;
	g0 = qp->p;
	xl = (double *) malloc(sizeof(double)*n);
	for (i=0;i<n;i++)
		xl[i] = 0;

	fatol = 0;
	frtol = 1e-12;
	fmin = -1e+32;
	cgtol = 0.1;	
	gtol = qp->eps;
	
	dtron(n, x, xl, xu, gtol, frtol, fatol, fmin, maxfev, cgtol); 

	free(xl);
}
