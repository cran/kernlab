extern double mymin(double, double);
extern double mymax(double, double); 

void dbreakpt(int n, double *x, double *xl, double *xu, double *w, int *nbrpt, double *brptmin, double *brptmax)
{
/*
c     **********
c
c     Subroutine dbreakpt
c
c     This subroutine computes the number of break-points, and
c     the minimal and maximal break-points of the projection of
c     x + alpha*w on the n-dimensional interval [xl,xu].
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
c       w is a double precision array of dimension n.
c         On entry w specifies the vector w.
c         On exit w is unchanged.
c
c       nbrpt is an integer variable.
c         On entry nbrpt need not be specified.
c         On exit nbrpt is the number of break points.
c
c       brptmin is a double precision variable
c         On entry brptmin need not be specified.
c         On exit brptmin is minimal break-point.
c
c       brptmax is a double precision variable
c         On entry brptmax need not be specified.
c         On exit brptmax is maximal break-point.
c
c     **********
*/
	int i;
	double brpt;
	
	*nbrpt = 0;
	for (i=0;i<n;i++)
		if (x[i] < xu[i] && w[i] > 0)
		{
			(*nbrpt)++;
			brpt = (xu[i] - x[i])/w[i];
			if (*nbrpt == 1)
				*brptmin = *brptmax = brpt;
			else
			{
				*brptmin = mymin(brpt, *brptmin);
				*brptmax = mymax(brpt, *brptmax);
			}
		}
		else
			if (x[i] > xl[i] && w[i] < 0)
			{
				(*nbrpt)++;
				brpt = (xl[i] - x[i])/w[i];
				if (*nbrpt == 1)
					*brptmin = *brptmax = brpt;
				else
				{
					*brptmin = mymin(brpt, *brptmin);
					*brptmax = mymax(brpt, *brptmax);
				}
			}
	if (*nbrpt == 0)
		*brptmin = *brptmax = 0;
}
