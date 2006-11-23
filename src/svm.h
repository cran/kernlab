#ifndef _LIBSVM_H
#define _LIBSVM_H

#ifdef __cplusplus
extern "C" {
#endif

struct svm_node
{
	int index;
	double value;
};

struct svm_problem
{
  int l, n;
	double *y;
	struct svm_node **x;
};

enum { C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR, C_BSVC, EPSILON_BSVR, SPOC, KBB };	/* svm_type */
enum { LINEAR, POLY, RBF, SIGMOID, R, LAPLACE, BESSEL, ANOVA, SPLINE };	/* kernel_type */

struct svm_parameter
{
	int svm_type;
	int kernel_type;
	int degree;	/* for poly */
	double gamma;	/* for poly/rbf/sigmoid */
	double coef0;	/* for poly/sigmoid */

	/* these are for training only */
	double cache_size; /* in MB */
	double eps;	/* stopping criteria */
	double C;	/* for C_SVC, EPSILON_SVR and NU_SVR */
	int nr_weight;		/* for C_SVC */
	int *weight_label;	/* for C_SVC */
	double* weight;		/* for C_SVC */
	double nu;	/* for NU_SVC, ONE_CLASS, and NU_SVR */
	double p;	/* for EPSILON_SVR */
	int shrinking;	/* use the shrinking heuristics */
        int qpsize;
        double Cbegin, Cstep;   /* for linear kernel */
        double lim; /* for bessel kernel */
        double *K; /* pointer to kernel matrix */
        int m;
};

struct BQP
 {
   double eps;
   int n;
   double *x, *C, *Q, *p;
};


#ifdef __cplusplus
}
#endif

#endif /* _LIBSVM_H */
