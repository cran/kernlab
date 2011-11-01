#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <R_ext/Lapack.h>
/* LAPACK */
/* extern int dpotf2_(char *, int *, double *, int *, int *); */

double dcholfact(int n, double *A, double *L)
{
	/* if A is p.d. , A = L*L'
	   if A is p.s.d. , A + lambda*I = L*L'; */  
	int indef, i;
	static double lambda = 1e-3/512/512;
	memcpy(L, A, sizeof(double)*n*n);
	F77_CALL(dpotf2)("L", &n, L, &n, &indef);
	if (indef != 0)
	{
		memcpy(L, A, sizeof(double)*n*n);
		for (i=0;i<n;i++)
			L[i*n+i] += lambda; 
		F77_CALL(dpotf2)("L", &n, L, &n, &indef);
		if (indef != 0)
		{
			//printf("A is not positive semi-definite\n");
			lambda *= 2;
		}
		return lambda;
	}
	return 0;
}

double dprecond(int n, double *A, double *C)
{
	/* Given a dense symmetric positive semidefinite matrix A, this
	subroutine computes the precondictioner C. Use full Cholesky
	factorization instead of incomplete Cholesky factorization that
	orginal TRON uses. It is the major difference between the two. */
	return dcholfact(n, A, C);
}
