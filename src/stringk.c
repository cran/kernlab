#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>


double ***cache ;


double kaux (const char *u, int p, const char *v, int q, int n, double lambda) {
  register int         j;
  double               tmp;

  /* case 1: if a full substring length is processed, return*/
  if (n == 0) return (1.0);

  /* check, if the value was already computed in a previous computation */
  if (cache [n] [p] [q] != -1.0) return (cache [n] [p] [q]); 
  
  /* case 2: at least one substring is to short */
  if (p < n || q < n) return (0.0);
    
  /* case 3: recursion */
  for (j= 0, tmp = 0; j < q; j++) {
    if (v [j] == u [p - 1]) 
      tmp += kaux (u, p - 1, v, j, n - 1, lambda) *   
        pow (lambda, (float) (q - j + 1));
  }

  cache [n] [p] [q] = lambda * kaux (u, p - 1, v, q, n, lambda) + tmp;
  return (cache [n] [p] [q]);
}


double seqk (const char *u, int p, const char *v, int q, int n, double lambda) {
  register int  j;
  double        kp;

  /* the simple case: (at least) one string is to short */
  if (p < n || q < n) return (0.0);

  /* the recursion: use kaux for the t'th substrings*/
  for (j = 0, kp = 0.0; j < q; j++) {
    if (v [j] == u [p - 1]) 
      kp += kaux (u, p - 1, v, j, n - 1, lambda) * lambda * lambda;
  }
  
  return (seqk (u, p - 1, v, q, n, lambda) + kp);
}

/* recursively computes the subsequence kernel between s1 and s2
   where subsequences of exactly length n are considered */

SEXP subsequencek(SEXP s1, SEXP s2, SEXP l1, SEXP l2, SEXP nr, SEXP lambdar) {

  const char *u = CHAR(STRING_ELT(s1, 0)); 
  const char *v = CHAR(STRING_ELT(s2, 0)); 
  int   p = *INTEGER(l1);
  int   q = *INTEGER(l2);
  int   n = *INTEGER(nr);
  double lambda = *REAL(lambdar);
  int           i, j, k;
  SEXP        ret;

  /* allocate memory for auxiallary cache variable */
  cache  = (double ***) malloc (n * sizeof (double **));
  for (i = 1; i < n; i++) {
    cache  [i] = (double **) malloc (p * sizeof (double *));
    for (j = 0; j < p; j++) {
      cache  [i] [j] = (double *) malloc (q * sizeof (double));
      for (k = 0; k < q; k++) 
	cache  [i] [j] [k] = -1.0;
    }
  }
  
  PROTECT(ret = allocVector(REALSXP, 1));

  /* invoke recursion */
  REAL(ret)[0] =  seqk (u, p, v, q, n, lambda);
 
  /* free memory */
  for (i = 1; i < n; i++) {
    for (j = 0; j < p; j++) 
      free (cache  [i] [j]);
    free (cache  [i]);
  }
  free (cache);
  UNPROTECT(1);

  return (ret);
}




/* computes the substring kernel between s1 and s2
   where substrings up to length n are considered */

SEXP fullsubstringk (SEXP s1, SEXP s2, SEXP l1, SEXP l2, SEXP nr, SEXP lambdar) {

  const char *u = CHAR(STRING_ELT(s1, 0)); 
  const char *v = CHAR(STRING_ELT(s2, 0)); 
  int   p = *INTEGER(l1);
  int   q = *INTEGER(l2);
  int   n = *INTEGER(nr);
  double lambda = *REAL(lambdar);
  register int  i, j, k;
  double        ret, tmp;
  SEXP          retk;

  /* computes the substring kernel */
  for (ret = 0.0, i = 0; i < p; i++) {
    for (j = 0; j < q; j++) 
      if (u [i] == v [j]) {
	for (k = 0, tmp = lambda * lambda;     /* starting condition */
	     (i + k < p) && (j + k < q) &&
	       (u [i + k] == v [j + k]) &&
	       (k < n);                        /* stop conditions */
	     k++, tmp *= (lambda * lambda))    /* update per iteration */
	  ret += tmp;
      }
  }
  
  PROTECT(retk = allocVector(REALSXP, 1));
  REAL(retk)[0] = ret; 
  UNPROTECT(1); 

  return (retk);
}

/* computes the substring kernel between s1 and s2
   where substrings of exactly length n are considered */

SEXP substringk (SEXP s1, SEXP s2, SEXP l1, SEXP l2, SEXP nr, SEXP lambdar) {

  const char *u = CHAR(STRING_ELT(s1, 0)); 
  const char *v = CHAR(STRING_ELT(s2, 0)); 
  int   p = *INTEGER(l1);
  int   q = *INTEGER(l2);
  int   n = *INTEGER(nr);
  double lambda = *REAL(lambdar);
  SEXP   retk;
 
  register int  i, j, k;
  double        ret, tmp;

  /* computes the substring kernel */
  for (ret = 0.0, i = 0; i < p; i++) {
    for (j = 0; j < q; j++) {
      for (k = 0, tmp = lambda * lambda;     /* starting condition */
	   (i + k < p) && (j + k < q) &&
	     (u [i + k] == v [j + k]) &&
	     (k < n);                        /* stop conditions */
	   k++, tmp *= (lambda * lambda));   /* update per iteration */
      
      if (k == n) ret += tmp;                /* update features in
						case of full match */
    }
  }
  
  PROTECT(retk = allocVector(REALSXP, 1));
  REAL(retk)[0] = ret; 
  UNPROTECT(1); 
  
  return (retk);
}
