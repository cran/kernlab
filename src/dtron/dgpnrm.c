#include <math.h>

double dgpnrm(int n, double *x, double *xl, double *xu, double *g)
{
/*
c     **********
c
c     Function dgpnrm
c
c     This function computes the infinite norm of the
c     projected gradient at x.
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
c       g is a double precision array of dimension n.
c         On entry g specifies the gradient g.
c         On exit g is unchanged.
c
c     **********
*/
	int i;
	double norm = 0;	

	for (i=0;i<n;i++)
		if (xl[i] != xu[i])
			if (!((g[i] <= 0 && x[i] == xu[i]) || (g[i] >= 0 && x[i] == xl[i])))
				if (fabs(g[i]) > norm)
					norm = fabs(g[i]);
	return norm;
}
