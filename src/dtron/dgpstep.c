void dgpstep(int n, double *x, double *xl, double *xu, double alpha, double *w, double *s)
{
/*
c     **********
c
c     Subroutine dgpstep
c
c     This subroutine computes the gradient projection step
c
c           s = P[x + alpha*w] - x,
c
c     where P is the projection on the n-dimensional interval [xl,xu].
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
c       alpha is a double precision variable.
c         On entry alpha specifies the scalar alpha.
c         On exit alpha is unchanged.
c
c       w is a double precision array of dimension n.
c         On entry w specifies the vector w.
c         On exit w is unchanged.
c
c       s is a double precision array of dimension n.
c         On entry s need not be specified.
c         On exit s contains the gradient projection step.
c
c     **********
*/
	int i;
	
	for (i=0;i<n;i++)
		if (x[i] + alpha*w[i] < xl[i])
			s[i] = xl[i] - x[i];
		else
			if (x[i] + alpha*w[i] > xu[i])
				s[i] = xu[i] - x[i];
			else
				s[i] = alpha*w[i];
}
