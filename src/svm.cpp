#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include "svm.h"
typedef float Qfloat;
typedef signed char schar;
#ifndef min
template <class T> inline T min(T x,T y) { return (x<y)?x:y; }
#endif
#ifndef max
template <class T> inline T max(T x,T y) { return (x>y)?x:y; }
#endif
template <class T> inline void swap(T& x, T& y) { T t=x; x=y; y=t; }
template <class S, class T> inline void clone(T*& dst, S* src, int n)
{
	dst = new T[n];
	memcpy((void *)dst,(void *)src,sizeof(T)*n);
}
#define INF DBL_MAX
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))


//

// Kernel Kcache
//
// l is the number of total data items
// size is the cache size limit in bytes
//
class Kcache
{
public:
	Kcache(int l,int size);
	~Kcache();

	// request data [0,len)
	// return some position p where [p,len) need to be filled
	// (p >= len if nothing needs to be filled)
	int get_data(const int index, Qfloat **data, int len);
	void swap_index(int i, int j);	// future_option
private:
	int l;
	int size;
	struct head_t
	{
		head_t *prev, *next;	// a cicular list
		Qfloat *data;
		int len;		// data[0,len) is cached in this entry
	};

	head_t* head;
	head_t lru_head;
	void lru_delete(head_t *h);
	void lru_insert(head_t *h);
};

Kcache::Kcache(int l_,int size_):l(l_),size(size_)
{
	head = (head_t *)calloc(l,sizeof(head_t));	// initialized to 0
	size /= sizeof(Qfloat);
	size -= l * sizeof(head_t) / sizeof(Qfloat);
	lru_head.next = lru_head.prev = &lru_head;
}

Kcache::~Kcache()
{
	for(head_t *h = lru_head.next; h != &lru_head; h=h->next)
		free(h->data);
	free(head);
}

void Kcache::lru_delete(head_t *h)
{
	// delete from current location
	h->prev->next = h->next;
	h->next->prev = h->prev;
}

void Kcache::lru_insert(head_t *h)
{
	// insert to last position
	h->next = &lru_head;
	h->prev = lru_head.prev;
	h->prev->next = h;
	h->next->prev = h;
}

int Kcache::get_data(const int index, Qfloat **data, int len)
{
	head_t *h = &head[index];
	if(h->len) lru_delete(h);
	int more = len - h->len;

	if(more > 0)
	{
		// free old space
		while(size < more)
		{ 
			head_t *old = lru_head.next;
			lru_delete(old);
			free(old->data);
			size += old->len;
			old->data = 0;
			old->len = 0;
		}

		// allocate new space
		h->data = (Qfloat *)realloc(h->data,sizeof(Qfloat)*len);
		size -= more;
		swap(h->len,len);
	}

	lru_insert(h);
	*data = h->data;
	return len;
}

void Kcache::swap_index(int i, int j)
{
	if(i==j) return;

	if(head[i].len) lru_delete(&head[i]);
	if(head[j].len) lru_delete(&head[j]);
	swap(head[i].data,head[j].data);
	swap(head[i].len,head[j].len);
	if(head[i].len) lru_insert(&head[i]);
	if(head[j].len) lru_insert(&head[j]);

	if(i>j) swap(i,j);
	for(head_t *h = lru_head.next; h!=&lru_head; h=h->next)
	{
		if(h->len > i)
		{
			if(h->len > j)
				swap(h->data[i],h->data[j]);
			else
			{
				// give up
				lru_delete(h);
				free(h->data);
				size += h->len;
				h->data = 0;
				h->len = 0;
			}
		}
	}
}

//
//Kernel evaluation
//
//the static method k_function is for doing single kernel evaluation
//the constructor of Kernel prepares to calculate the l*l kernel matrix
//the member function get_Q is for getting one column from the Q Matrix
//
class kernel {
public:
	kernel(int l, svm_node * const * x, const svm_parameter& param);
	virtual ~kernel();

	static double k_function(const svm_node *x, const svm_node *y,
				 const svm_parameter& param);
	virtual Qfloat *get_Q(int column, int len) const = 0;
	virtual void swap_index(int i, int j) const	// no so const...
	{
		swap(x[i],x[j]);
		if(x_square) swap(x_square[i],x_square[j]);
	}
protected:

	double (kernel::*kernel_function)(int i, int j) const;

private:
	const svm_node **x;
	double *x_square;
     
	// svm_parameter
	const int kernel_type;
	const double degree;
	const double gamma;
        const double coef0;
        const double lim;
        const int d;
      
                

	static double dot(const svm_node *px, const svm_node *py);
        static double anova(const svm_node *px, const svm_node *py, const double sigma, const int d, const double degree);

	double kernel_linear(int i, int j) const
	{
		return dot(x[i],x[j]);
	}
	double kernel_poly(int i, int j) const
	{
		return pow(gamma*dot(x[i],x[j])+coef0,degree);
	}
	double kernel_rbf(int i, int j) const
        {
         return exp(-gamma*(x_square[i]+x_square[j]-2*dot(x[i],x[j])));
        }
	double kernel_laplace(int i, int j) const
	{
	return exp(-gamma*sqrt(fabs(x_square[i]+x_square[j]-2*dot(x[i],x[j]))));
	}
        double kernel_bessel(int i, int j) const
        { 
	  double bkt = gamma*sqrt(fabs(x_square[i]+x_square[j]-2*dot(x[i],x[j])));
	  if (bkt < 0.000001){
            return 1 ;
	  }
          else {
	    return(pow(((jn((int)degree, bkt)*pow(bkt,(-degree)))/lim),coef0));
	  }
	}
	double kernel_sigmoid(int i, int j) const
	{
		return tanh(gamma*dot(x[i],x[j])+coef0);
	}
        double kernel_anova(int i, int j) const
        {
	 return  anova(x[i], x[j], gamma, d, degree);
	}
};

kernel::kernel(int l, svm_node * const * x_, const svm_parameter& param)
:kernel_type(param.kernel_type), degree(param.degree),
 gamma(param.gamma), coef0(param.coef0), lim(param.lim), d(param.d)
{
	switch(kernel_type)
	{
		case LINEAR:
			kernel_function = &kernel::kernel_linear;
			break;
		case POLY:
			kernel_function = &kernel::kernel_poly;
			break;
		case RBF:
			kernel_function = &kernel::kernel_rbf;
			break;
		case LAPLACE:
			kernel_function = &kernel::kernel_laplace;
			break;
        	case BESSEL:
	                kernel_function = &kernel::kernel_bessel;
			break;
		case SIGMOID:
			kernel_function = &kernel::kernel_sigmoid;
			break;
		case ANOVA:
		        kernel_function = &kernel::kernel_anova;
			break;

	}

	clone(x,x_,l);

	if(kernel_type == RBF || kernel_type == LAPLACE || kernel_type == BESSEL)
	{
		x_square = new double[l];
		for(int i=0;i<l;i++)
			x_square[i] = dot(x[i],x[i]);
	}

	else
		x_square = 0;
}

kernel::~kernel()
  {
	delete[] x;
	delete[] x_square;
  }

double kernel::dot(const svm_node *px, const svm_node *py)
{
	double sum = 0;
	while(px->index != -1 && py->index != -1)
	{
		if(px->index == py->index)
		  {
			sum += px->value * py->value;
			++px;
			++py;
		}
		else
		  { 
			if(px->index > py->index)
				++py;
			else
				++px;
		}			
	}
	return sum;
}

double kernel::anova(const svm_node *px, const svm_node *py, const double sigma, const int d, const double degree)
{

	double sum = 0;
        //int zero = d; 
	double tv;
	while(px->index != -1 && py->index != -1)
	{ 
		if(px->index == py->index)
		  { 	
		    tv = (px->value - py->value) * (px->value - py->value);
		    // if (fabs(tv) < 0.000001)
		    //	sum ++;
		    // else
			sum += exp( - sigma * tv);
		    ++px;
		    ++py;
		}
		else
		{
			if(px->index > py->index)
			  { 
			    sum += exp( - sigma * (py->value * py->value));
			    ++py;
			  }
			else
			  {
			    sum += exp( - sigma * (px->value * px->value));
			    ++px;
			  }
		}
	//	zero--;
	}
	return (pow(sum,degree));
}




double kernel::k_function(const svm_node *x, const svm_node *y,
			  const svm_parameter& param)
{
	switch(param.kernel_type)
	{
		case LINEAR:
			return dot(x,y);
		case POLY:
			return pow(param.gamma*dot(x,y)+param.coef0,param.degree);
		case RBF:
		{
			double sum = 0;
			while(x->index != -1 && y->index !=-1)
			{
				if(x->index == y->index)
				{
					double d = x->value - y->value;
					sum += d*d;
					++x;
					++y;
				}
				else
				{
					if(x->index > y->index)
					{	
						sum += y->value * y->value;
						++y;
					}
					else
					{
						sum += x->value * x->value;
						++x;
					}
				}
			}

			while(x->index != -1)
			{
				sum += x->value * x->value;
				++x;
			}

			while(y->index != -1)
			{
				sum += y->value * y->value;
				++y;
			}
			
			return exp(-param.gamma*sum);
		}
		case SIGMOID:
			return tanh(param.gamma*dot(x,y)+param.coef0);
		default:
			return 0;	/* Unreachable */
	}
}

//Generalized SMO+SVMlight algorithm
//Solves:
//
//	min 0.5(\alpha^T Q \alpha) + b^T \alpha
//
//		y^T \alpha = \delta
//		y_i = +1 or -1
//		0 <= alpha_i <= Cp for y_i = 1
//		0 <= alpha_i <= Cn for y_i = -1
//
//Given:
//
//	Q, b, y, Cp, Cn, and an initial feasible point \alpha
//	l is the size of vectors and matrices
//	eps is the stopping criterion
//
//solution will be put in \alpha, objective value will be put in obj
//
class Solver {
public:
	Solver() {};
	virtual ~Solver() {};

	struct SolutionInfo {
		double obj;
		double rho;
		double upper_bound_p;
		double upper_bound_n;
		double r;	// for Solver_NU
	};

	void Solve(int l, const kernel& Q, const double *b_, const schar *y_,
		   double *alpha_, double Cp, double Cn, double eps,
		   SolutionInfo* si, int shrinking);
protected:
	int active_size;
	schar *y;
	double *G;		// gradient of objective function
	enum { LOWER_BOUND, UPPER_BOUND, FREE };
	char *alpha_status;	// LOWER_BOUND, UPPER_BOUND, FREE
	double *alpha;
	const kernel *Q;
	double eps;
	double Cp,Cn;
	double *b;
	int *active_set;
	double *G_bar;		// gradient, if we treat free variables as 0
	int l;
	bool unshrinked;	// XXX

	double get_C(int i)
	{
		return (y[i] > 0)? Cp : Cn;
	}
	void update_alpha_status(int i)
	{
		if(alpha[i] >= get_C(i))
			alpha_status[i] = UPPER_BOUND;
		else if(alpha[i] <= 0)
			alpha_status[i] = LOWER_BOUND;
		else alpha_status[i] = FREE;
	}
	bool is_upper_bound(int i) { return alpha_status[i] == UPPER_BOUND; }
	bool is_lower_bound(int i) { return alpha_status[i] == LOWER_BOUND; }
	bool is_free(int i) { return alpha_status[i] == FREE; }
	void swap_index(int i, int j);
	void reconstruct_gradient();
	virtual int select_working_set(int &i, int &j);
	virtual double calculate_rho();
	virtual void do_shrinking();
};

void Solver::swap_index(int i, int j)
{
	Q->swap_index(i,j);
	swap(y[i],y[j]);
	swap(G[i],G[j]);
	swap(alpha_status[i],alpha_status[j]);
	swap(alpha[i],alpha[j]);
	swap(b[i],b[j]);
	swap(active_set[i],active_set[j]);
	swap(G_bar[i],G_bar[j]);
}

void Solver::reconstruct_gradient()
{
	// reconstruct inactive elements of G from G_bar and free variables

	if(active_size == l) return;

	int i;
	for(i=active_size;i<l;i++)
		G[i] = G_bar[i] + b[i];
	
	for(i=0;i<active_size;i++)
		if(is_free(i))
		{
			const Qfloat *Q_i = Q->get_Q(i,l);
			double alpha_i = alpha[i];
			for(int j=active_size;j<l;j++)
				G[j] += alpha_i * Q_i[j];
		}
}

void Solver::Solve(int l, const kernel& Q, const double *b_, const schar *y_,
		   double *alpha_, double Cp, double Cn, double eps,
		   SolutionInfo* si, int shrinking)
{
	this->l = l;
	this->Q = &Q;
	clone(b, b_,l);
	clone(y, y_,l);
	clone(alpha,alpha_,l);
	this->Cp = Cp;
	this->Cn = Cn;
	this->eps = eps;
	unshrinked = false;

	// initialize alpha_status
	{
		alpha_status = new char[l];
		for(int i=0;i<l;i++)
			update_alpha_status(i);
	}

	// initialize active set (for shrinking)
	{
		active_set = new int[l];
		for(int i=0;i<l;i++)
			active_set[i] = i;
		active_size = l;
	}

	// initialize gradient
	{
		G = new double[l];
		G_bar = new double[l];
		int i;
		for(i=0;i<l;i++)
		{
			G[i] = b[i];
			G_bar[i] = 0;
		}
		for(i=0;i<l;i++)
			if(!is_lower_bound(i))
			{
				Qfloat *Q_i = Q.get_Q(i,l);
				double alpha_i = alpha[i];
				int j;
				for(j=0;j<l;j++)
					G[j] += alpha_i*Q_i[j];
				if(is_upper_bound(i))
					for(j=0;j<l;j++)
						G_bar[j] += get_C(i) * Q_i[j];
			}
	}

	// optimization step

	int iter = 0;
	int counter = min(l,1000)+1;

	while(1)
	{
		// show progress and do shrinking
	 
		if(--counter == 0)
		{ 
			counter = min(l,1000);
			if(shrinking) do_shrinking();
			//	info("."); info_flush();
		}

		int i,j;
		if(select_working_set(i,j)!=0)
		{
			// reconstruct the whole gradient
			reconstruct_gradient();
			// reset active set size and check
			active_size = l;
			//info("*"); info_flush();
			if(select_working_set(i,j)!=0)
				break;
			else
				counter = 1;	// do shrinking next iteration
		}
		
		++iter;

		// update alpha[i] and alpha[j], handle bounds carefully
		
		const Qfloat *Q_i = Q.get_Q(i,active_size);
		const Qfloat *Q_j = Q.get_Q(j,active_size);

		double C_i = get_C(i);
		double C_j = get_C(j);

		double old_alpha_i = alpha[i];
		double old_alpha_j = alpha[j];

		if(y[i]!=y[j])
		{
			double delta = (-G[i]-G[j])/max(Q_i[i]+Q_j[j]+2*Q_i[j],(Qfloat)0);
			double diff = alpha[i] - alpha[j];
			alpha[i] += delta;
			alpha[j] += delta;
			
			if(diff > 0)
			{
				if(alpha[j] < 0)
				{
					alpha[j] = 0;
					alpha[i] = diff;
				}
			}
			else
			{
				if(alpha[i] < 0)
				{
					alpha[i] = 0;
					alpha[j] = -diff;
				}
			}
			if(diff > C_i - C_j)
			{
				if(alpha[i] > C_i)
				{
					alpha[i] = C_i;
					alpha[j] = C_i - diff;
				}
			}
			else
			{
				if(alpha[j] > C_j)
				{
					alpha[j] = C_j;
					alpha[i] = C_j + diff;
				}
			}
		}
		else
		{
			double delta = (G[i]-G[j])/max(Q_i[i]+Q_j[j]-2*Q_i[j],(Qfloat)0);
			double sum = alpha[i] + alpha[j];
			alpha[i] -= delta;
			alpha[j] += delta;
			if(sum > C_i)
			{
				if(alpha[i] > C_i)
				{
					alpha[i] = C_i;
					alpha[j] = sum - C_i;
				}
			}
			else
			{
				if(alpha[j] < 0)
				{
					alpha[j] = 0;
					alpha[i] = sum;
				}
			}
			if(sum > C_j)
			{
				if(alpha[j] > C_j)
				{
					alpha[j] = C_j;
					alpha[i] = sum - C_j;
				}
			}
			else
			{
				if(alpha[i] < 0)
				{
					alpha[i] = 0;
					alpha[j] = sum;
				}
			}
		}

		// update G

		double delta_alpha_i = alpha[i] - old_alpha_i;
		double delta_alpha_j = alpha[j] - old_alpha_j;
		
		for(int k=0;k<active_size;k++)
		{
			G[k] += Q_i[k]*delta_alpha_i + Q_j[k]*delta_alpha_j;
		}

		// update alpha_status and G_bar

		{
			bool ui = is_upper_bound(i);
			bool uj = is_upper_bound(j);
			update_alpha_status(i);
			update_alpha_status(j);
			int k;
			if(ui != is_upper_bound(i))
			{
				Q_i = Q.get_Q(i,l);
				if(ui)
					for(k=0;k<l;k++)
						G_bar[k] -= C_i * Q_i[k];
				else
					for(k=0;k<l;k++)
						G_bar[k] += C_i * Q_i[k];
			}

			if(uj != is_upper_bound(j))
			{
				Q_j = Q.get_Q(j,l);
				if(uj)
					for(k=0;k<l;k++)
						G_bar[k] -= C_j * Q_j[k];
				else
					for(k=0;k<l;k++)
						G_bar[k] += C_j * Q_j[k];
			}
		}
	}

	// calculate rho

	si->rho = calculate_rho();

	// calculate objective value
	{
		double v = 0;
		int i;
		for(i=0;i<l;i++)
			v += alpha[i] * (G[i] + b[i]);

		si->obj = v/2;
	}

	// put back the solution
	{
		for(int i=0;i<l;i++)
			alpha_[active_set[i]] = alpha[i];
	}

	// juggle everything back
	/*{
		for(int i=0;i<l;i++)
			while(active_set[i] != i)
				swap_index(i,active_set[i]);
				// or Q.swap_index(i,active_set[i]);
	}*/

	si->upper_bound_p = Cp;
	si->upper_bound_n = Cn;

	//info("\noptimization finished, #iter = %d\n",iter);

	delete[] b;
	delete[] y;
	delete[] alpha;
	delete[] alpha_status;
	delete[] active_set;
	delete[] G;
	delete[] G_bar;
}

// return 1 if already optimal, return 0 otherwise
int Solver::select_working_set(int &out_i, int &out_j)
{
	// return i,j which maximize -grad(f)^T d , under constraint
	// if alpha_i == C, d != +1
	// if alpha_i == 0, d != -1

	double Gmax1 = -INF;		// max { -grad(f)_i * d | y_i*d = +1 }
	int Gmax1_idx = -1;

	double Gmax2 = -INF;		// max { -grad(f)_i * d | y_i*d = -1 }
	int Gmax2_idx = -1;

	for(int i=0;i<active_size;i++)
	{
		if(y[i]==+1)	// y = +1
		{
			if(!is_upper_bound(i))	// d = +1
			{
				if(-G[i] > Gmax1)
				{
					Gmax1 = -G[i];
					Gmax1_idx = i;
				}
			}
			if(!is_lower_bound(i))	// d = -1
			{
				if(G[i] > Gmax2)
				{
					Gmax2 = G[i];
					Gmax2_idx = i;
				}
			}
		}
		else		// y = -1
		{
			if(!is_upper_bound(i))	// d = +1
			{
				if(-G[i] > Gmax2)
				{
					Gmax2 = -G[i];
					Gmax2_idx = i;
				}
			}
			if(!is_lower_bound(i))	// d = -1
			{
				if(G[i] > Gmax1)
				{
					Gmax1 = G[i];
					Gmax1_idx = i;
				}
			}
		}
	}
	if(Gmax1+Gmax2 < eps)
		return 1;

	out_i = Gmax1_idx;
	out_j = Gmax2_idx;
	return 0;
}

void Solver::do_shrinking()
{
	int i,j,k;
	if(select_working_set(i,j)!=0) return;
	double Gm1 = -y[j]*G[j];
	double Gm2 = y[i]*G[i];

	// shrink
	
	for(k=0;k<active_size;k++)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] >= Gm1) continue;
			}
			else	if(-G[k] >= Gm2) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] >= Gm2) continue;
			}
			else	if(G[k] >= Gm1) continue;
		}
		else continue;

		--active_size;
		swap_index(k,active_size);
		--k;	// look at the newcomer
	}

	// unshrink, check all variables again before final iterations

	if(unshrinked || -(Gm1 + Gm2) > eps*10) return;
	
	unshrinked = true;
	reconstruct_gradient();

	for(k=l-1;k>=active_size;k--)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] < Gm1) continue;
			}
			else	if(-G[k] < Gm2) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] < Gm2) continue;
			}
			else	if(G[k] < Gm1) continue;
		}
		else continue;

		swap_index(k,active_size);
		active_size++;
		++k;	// look at the newcomer
	}
}

double Solver::calculate_rho()
{
	double r;
	int nr_free = 0;
	double ub = INF, lb = -INF, sum_free = 0;
	for(int i=0;i<active_size;i++)
	{
		double yG = y[i]*G[i];

		if(is_lower_bound(i))
		{
			if(y[i] > 0)
				ub = min(ub,yG);
			else
				lb = max(lb,yG);
		}
		else if(is_upper_bound(i))
		{
			if(y[i] < 0)
				ub = min(ub,yG);
			else
				lb = max(lb,yG);
		}
		else
		{
			++nr_free;
			sum_free += yG;
		}
	}

	if(nr_free>0)
		r = sum_free/nr_free;
	else
		r = (ub+lb)/2;

	return r;
}

//
// Solver for nu-svm classification and regression
//
// additional constraint: e^T \alpha = constant
//
class Solver_NU : public Solver
{
public:
	Solver_NU() {}
	void Solve(int l, const kernel& Q, const double *b, const schar *y,
		   double *alpha, double Cp, double Cn, double eps,
		   SolutionInfo* si, int shrinking)
	{
		this->si = si;
		Solver::Solve(l,Q,b,y,alpha,Cp,Cn,eps,si,shrinking);
	}
private:
	SolutionInfo *si;
	int select_working_set(int &i, int &j);
	double calculate_rho();
	void do_shrinking();
};

int Solver_NU::select_working_set(int &out_i, int &out_j)
{
	// return i,j which maximize -grad(f)^T d , under constraint
	// if alpha_i == C, d != +1
	// if alpha_i == 0, d != -1

	double Gmax1 = -INF;	// max { -grad(f)_i * d | y_i = +1, d = +1 }
	int Gmax1_idx = -1;

	double Gmax2 = -INF;	// max { -grad(f)_i * d | y_i = +1, d = -1 }
	int Gmax2_idx = -1;

	double Gmax3 = -INF;	// max { -grad(f)_i * d | y_i = -1, d = +1 }
	int Gmax3_idx = -1;

	double Gmax4 = -INF;	// max { -grad(f)_i * d | y_i = -1, d = -1 }
	int Gmax4_idx = -1;

	for(int i=0;i<active_size;i++)
	{
		if(y[i]==+1)	// y == +1
		{
			if(!is_upper_bound(i))	// d = +1
			{
				if(-G[i] > Gmax1)
				{
					Gmax1 = -G[i];
					Gmax1_idx = i;
				}
			}
			if(!is_lower_bound(i))	// d = -1
			{
				if(G[i] > Gmax2)
				{
					Gmax2 = G[i];
					Gmax2_idx = i;
				}
			}
		}
		else		// y == -1
		{
			if(!is_upper_bound(i))	// d = +1
			{
				if(-G[i] > Gmax3)
				{
					Gmax3 = -G[i];
					Gmax3_idx = i;
				}
			}
			if(!is_lower_bound(i))	// d = -1
			{
				if(G[i] > Gmax4)
				{
					Gmax4 = G[i];
					Gmax4_idx = i;
				}
			}
		}
	}

	if(max(Gmax1+Gmax2,Gmax3+Gmax4) < eps)
 		return 1;

	if(Gmax1+Gmax2 > Gmax3+Gmax4)
	{
		out_i = Gmax1_idx;
		out_j = Gmax2_idx;
	}
	else
	{
		out_i = Gmax3_idx;
		out_j = Gmax4_idx;
	}
	return 0;
}

void Solver_NU::do_shrinking()
{
	double Gmax1 = -INF;	// max { -grad(f)_i * d | y_i = +1, d = +1 }
	double Gmax2 = -INF;	// max { -grad(f)_i * d | y_i = +1, d = -1 }
	double Gmax3 = -INF;	// max { -grad(f)_i * d | y_i = -1, d = +1 }
	double Gmax4 = -INF;	// max { -grad(f)_i * d | y_i = -1, d = -1 }

	int k;
	for(k=0;k<active_size;k++)
	{
		if(!is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] > Gmax1) Gmax1 = -G[k];
			}
			else	if(-G[k] > Gmax3) Gmax3 = -G[k];
		}
		if(!is_lower_bound(k))
		{
			if(y[k]==+1)
			{	
				if(G[k] > Gmax2) Gmax2 = G[k];
			}
			else	if(G[k] > Gmax4) Gmax4 = G[k];
		}
	}

	double Gm1 = -Gmax2;
	double Gm2 = -Gmax1;
	double Gm3 = -Gmax4;
	double Gm4 = -Gmax3;

	for(k=0;k<active_size;k++)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] >= Gm1) continue;
			}
			else	if(-G[k] >= Gm3) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] >= Gm2) continue;
			}
			else	if(G[k] >= Gm4) continue;
		}
		else continue;

		--active_size;
		swap_index(k,active_size);
		--k;	// look at the newcomer
	}

	// unshrink, check all variables again before final iterations

	if(unshrinked || max(-(Gm1+Gm2),-(Gm3+Gm4)) > eps*10) return;
	
	unshrinked = true;
	reconstruct_gradient();

	for(k=l-1;k>=active_size;k--)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] < Gm1) continue;
			}
			else	if(-G[k] < Gm3) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] < Gm2) continue;
			}
			else	if(G[k] < Gm4) continue;
		}
		else continue;

		swap_index(k,active_size);
		active_size++;
		++k;	// look at the newcomer
	}
}

double Solver_NU::calculate_rho()
{
	int nr_free1 = 0,nr_free2 = 0;
	double ub1 = INF, ub2 = INF;
	double lb1 = -INF, lb2 = -INF;
	double sum_free1 = 0, sum_free2 = 0;

	for(int i=0;i<active_size;i++)
	{
		if(y[i]==+1)
		{
			if(is_lower_bound(i))
				ub1 = min(ub1,G[i]);
			else if(is_upper_bound(i))
				lb1 = max(lb1,G[i]);
			else
			{
				++nr_free1;
				sum_free1 += G[i];
			}
		}
		else
		{
			if(is_lower_bound(i))
				ub2 = min(ub2,G[i]);
			else if(is_upper_bound(i))
				lb2 = max(lb2,G[i]);
			else
			{
				++nr_free2;
				sum_free2 += G[i];
			}
		}
	}

	double r1,r2;
	if(nr_free1 > 0)
		r1 = sum_free1/nr_free1;
	else
		r1 = (ub1+lb1)/2;
	
	if(nr_free2 > 0)
		r2 = sum_free2/nr_free2;
	else
		r2 = (ub2+lb2)/2;
	
	si->r = (r1+r2)/2;
	return (r1-r2)/2;
}

//
// Q matrices for various formulations
//
class SVC_Q: public kernel
{ 
public:
	SVC_Q(const svm_problem& prob, const svm_parameter& param, const schar *y_)
	:kernel(prob.l, prob.x, param)
	{
		clone(y,y_,prob.l);
		cache = new Kcache(prob.l,(int)(param.cache_size*(1<<20)));
		ktype = param.kernel_type;
		expr = param.expr;
		rho = param.rho;

	}
	
	Qfloat *get_Q(int i, int len) const
	{
		Qfloat *data;
		int start; 
	
		if((start = cache->get_data(i,&data,len) < len)  && ktype != 4)
		  { 
		    for(int j=start;j<len;j++)
		      data[j] = (Qfloat)(y[i]*y[j]*(this->*kernel_function)(i,j));
		  }
		else if(ktype == 4)
		  {
		   SEXP R_fcall, dummy, ans;
		   PROTECT(R_fcall = lang2(expr, R_NilValue));
		   PROTECT(ans= allocSExp(REALSXP));
		   PROTECT(dummy = allocVector(INTSXP, 1));
		   INTEGER(dummy)[0] = i+1;
		   SETCADR(R_fcall,dummy);
		    //	   SET_VECTOR_ELT(ans,0,eval(R_fcall,rho));
		    //	   ans=eval(R_fcall,rho);
		   ans = eval(R_fcall,rho);
		   for(int j=start;j<len;j++)
		     data[j] = (Qfloat)(y[i]*y[j]*REAL(ans)[j]);
		   UNPROTECT(3);
		}
		return data;
	}

	void swap_index(int i, int j) const
	{
		cache->swap_index(i,j);
		kernel::swap_index(i,j);
		swap(y[i],y[j]);
	}

	~SVC_Q()
	{
		delete[] y;
		delete cache;
	}
private:
        SEXP expr;
        SEXP rho;
        int ktype;
        schar *y;
	Kcache *cache;
};

class ONE_CLASS_Q: public kernel
{
public:
	ONE_CLASS_Q(const svm_problem& prob, const svm_parameter& param)
	:kernel(prob.l, prob.x, param)
	{
		cache = new Kcache(prob.l,(int)(param.cache_size*(1<<20)));
		ktype = param.kernel_type;
		expr = param.expr;
		rho = param.rho;
	}
	
	Qfloat *get_Q(int i, int len) const
	{
		Qfloat *data;
		int start;
		if((start = cache->get_data(i,&data,len)) < len && ktype !=4)
		{
			for(int j=start;j<len;j++)
				data[j] = (Qfloat)(this->*kernel_function)(i,j);
		}
			else if(ktype == 4)
		  {
		   SEXP R_fcall, dummy, ans;
		   PROTECT(R_fcall = lang2(expr, R_NilValue));
		   PROTECT(ans= allocSExp(REALSXP));
		   PROTECT(dummy = allocVector(INTSXP, 1));
		   INTEGER(dummy)[0] = i+1;
		   SETCADR(R_fcall,dummy);
		   ans = eval(R_fcall,rho);
		   for(int j=start;j<len;j++)
		     data[j] = (Qfloat)(REAL(ans)[j]);
		   UNPROTECT(3);
		  }
		return data;
	}

	void swap_index(int i, int j) const
	{
		cache->swap_index(i,j);
		kernel::swap_index(i,j);
	}

	~ONE_CLASS_Q()
	{
		delete cache;
	}
private:
	Kcache *cache;
        SEXP expr;
        SEXP rho;
        int ktype;
};

class SVR_Q: public kernel
{ 
public:
	SVR_Q(const svm_problem& prob, const svm_parameter& param)
	:kernel(prob.l, prob.x, param)
	{
		l = prob.l;
		cache = new Kcache(l,(int)(param.cache_size*(1<<20)));
		sign = new schar[2*l];
		index = new int[2*l];
		ktype = param.kernel_type;
		expr = param.expr;
		rho = param.rho;
		for(int k=0;k<l;k++)
		{
			sign[k] = 1;
			sign[k+l] = -1;
			index[k] = k;
			index[k+l] = k;
		}
		buffer[0] = new Qfloat[2*l];
		buffer[1] = new Qfloat[2*l];
		next_buffer = 0;
	}

	void swap_index(int i, int j) const
	{
		swap(sign[i],sign[j]);
		swap(index[i],index[j]);
	}
	
	Qfloat *get_Q(int i, int len) const
	{
		Qfloat *data;
		int real_i = index[i];
		if(cache->get_data(real_i,&data,l) < l)
		{
			for(int j=0;j<l;j++)
				data[j] = (Qfloat)(this->*kernel_function)(real_i,j);
		}
		else if(ktype == 4)
		  {
		    SEXP R_fcall, dummy, ans;
		    PROTECT(R_fcall = lang2(expr, R_NilValue));
		    PROTECT(ans= allocSExp(REALSXP));
		    PROTECT(dummy = allocVector(INTSXP, 1));
		    INTEGER(dummy)[0] = real_i+1;
		    SETCADR(R_fcall,dummy);
		    ans = eval(R_fcall,rho);
		    for(int j=0;j<l;j++)
		      data[j] = (Qfloat)(REAL(ans)[j]);
		    UNPROTECT(3);
		  }
	
		// reorder and copy
		Qfloat *buf = buffer[next_buffer];
		next_buffer = 1 - next_buffer;
		schar si = sign[i];
		for(int j=0;j<len;j++)
			buf[j] = si * sign[j] * data[index[j]];
		return buf;
	}

	~SVR_Q()
	{
		delete cache;
		delete[] sign;
		delete[] index;
		delete[] buffer[0];
		delete[] buffer[1];
	}
private:
	int l;
	Kcache *cache;
	schar *sign;
	int *index;
	mutable int next_buffer;
	Qfloat* buffer[2];
       SEXP expr;
        SEXP rho;
        int ktype;
};

//
// construct and solve various formulations
//
static void solve_c_svc(
	const svm_problem *prob, const svm_parameter* param,
	double *alpha, Solver::SolutionInfo* si, double Cp, double Cn)
{
	int l = prob->l;
	double *minus_ones = new double[l];
	schar *y = new schar[l];

	int i;

	for(i=0;i<l;i++)
	{
		alpha[i] = 0;
		minus_ones[i] = -1;
		if(prob->y[i] > 0) y[i] = +1; else y[i]=-1;
	}

	Solver s;
	s.Solve(l, SVC_Q(*prob,*param,y), minus_ones, y,
		alpha, Cp, Cn, param->eps, si, param->shrinking);

	double sum_alpha=0;
	for(i=0;i<l;i++)
		sum_alpha += alpha[i];

	//info("nu = %f\n", sum_alpha/(param->C*prob->l));

	for(i=0;i<l;i++)
		alpha[i] *= y[i];

	delete[] minus_ones;
	delete[] y;
}

static void solve_nu_svc(
	const svm_problem *prob, const svm_parameter *param,
	double *alpha, Solver::SolutionInfo* si)
{
	int i;
	int l = prob->l;
	double nu = param->nu;

	schar *y = new schar[l];

	for(i=0;i<l;i++)
		if(prob->y[i]>0)
			y[i] = +1;
		else
			y[i] = -1;

	double sum_pos = nu*l/2;
	double sum_neg = nu*l/2;

	for(i=0;i<l;i++)
		if(y[i] == +1)
		{
			alpha[i] = min(1.0,sum_pos);
			sum_pos -= alpha[i];
		}
		else
		{
			alpha[i] = min(1.0,sum_neg);
			sum_neg -= alpha[i];
		}

	double *zeros = new double[l];

	for(i=0;i<l;i++)
		zeros[i] = 0;

	Solver_NU s;
	s.Solve(l, SVC_Q(*prob,*param,y), zeros, y,
		alpha, 1.0, 1.0, param->eps, si,  param->shrinking);
	double r = si->r;

	//info("C = %f\n",1/r);

	for(i=0;i<l;i++)
		alpha[i] *= y[i]/r;

	si->rho /= r;
	si->obj /= (r*r);
	si->upper_bound_p = 1/r;
	si->upper_bound_n = 1/r;

	delete[] y;
	delete[] zeros;
}

static void solve_one_class(
	const svm_problem *prob, const svm_parameter *param,
	double *alpha, Solver::SolutionInfo* si)
{
	int l = prob->l;
	double *zeros = new double[l];
	schar *ones = new schar[l];
	int i;

	int n = (int)(param->nu*prob->l);	// # of alpha's at upper bound

	for(i=0;i<n;i++)
		alpha[i] = 1;
	alpha[n] = param->nu * prob->l - n;
	for(i=n+1;i<l;i++)
		alpha[i] = 0;

	for(i=0;i<l;i++)
	{
		zeros[i] = 0;
		ones[i] = 1;
	}

	Solver s;
	s.Solve(l, ONE_CLASS_Q(*prob,*param), zeros, ones,
		alpha, 1.0, 1.0, param->eps, si, param->shrinking);

	delete[] zeros;
	delete[] ones;
}

static void solve_epsilon_svr(
	const svm_problem *prob, const svm_parameter *param,
	double *alpha, Solver::SolutionInfo* si)
{
	int l = prob->l;
	double *alpha2 = new double[2*l];
	double *linear_term = new double[2*l];
	schar *y = new schar[2*l];
	int i;

	for(i=0;i<l;i++)
	{
		alpha2[i] = 0;
		linear_term[i] = param->p - prob->y[i];
		y[i] = 1;

		alpha2[i+l] = 0;
		linear_term[i+l] = param->p + prob->y[i];
		y[i+l] = -1;
	}

	Solver s;
	s.Solve(2*l, SVR_Q(*prob,*param), linear_term, y,
		alpha2, param->C, param->C, param->eps, si, param->shrinking);

	double sum_alpha = 0;
	for(i=0;i<l;i++)
	{
		alpha[i] = alpha2[i] - alpha2[i+l];
		sum_alpha += fabs(alpha[i]);
	}
	//info("nu = %f\n",sum_alpha/(param->C*l));

	delete[] alpha2;
	delete[] linear_term;
	delete[] y;
}

static void solve_nu_svr(
	const svm_problem *prob, const svm_parameter *param,
	double *alpha, Solver::SolutionInfo* si)
{
	int l = prob->l;
	double C = param->C;
	double *alpha2 = new double[2*l];
	double *linear_term = new double[2*l];
	schar *y = new schar[2*l];
	int i;

	double sum = C * param->nu * l / 2;
	for(i=0;i<l;i++)
	{
		alpha2[i] = alpha2[i+l] = min(sum,C);
		sum -= alpha2[i];

		linear_term[i] = - prob->y[i];
		y[i] = 1;

		linear_term[i+l] = prob->y[i];
		y[i+l] = -1;
	}

	Solver_NU s;
	s.Solve(2*l, SVR_Q(*prob,*param), linear_term, y,
		alpha2, C, C, param->eps, si, param->shrinking);

	//info("epsilon = %f\n",-si->r);

	for(i=0;i<l;i++)
		alpha[i] = alpha2[i] - alpha2[i+l];

	delete[] alpha2;
	delete[] linear_term;
	delete[] y;
}

//
// decision_function
//
struct decision_function
{
	double *alpha;
	double rho;	
};

decision_function svm_train_one(
	const svm_problem *prob, const svm_parameter *param,
	double Cp, double Cn)
{
	double *alpha = Malloc(double,prob->l);
	Solver::SolutionInfo si;
	switch(param->svm_type)
	{
		case C_SVC:
			solve_c_svc(prob,param,alpha,&si,Cp,Cn);
			break;
		case NU_SVC:
			solve_nu_svc(prob,param,alpha,&si);
			break;
		case ONE_CLASS:
			solve_one_class(prob,param,alpha,&si);
			break;
		case EPSILON_SVR:
			solve_epsilon_svr(prob,param,alpha,&si);
			break;
		case NU_SVR:
			solve_nu_svr(prob,param,alpha,&si);
			break;
	}

	//info("obj = %f, rho = %f\n",si.obj,si.rho);

	// output SVs

// 	int nSV = 0;
// 	int nBSV = 0;
// 	for(int i=0;i<prob->l;i++)
// 	{
// 		if(fabs(alpha[i]) > 0)
// 		{
// 			++nSV;
// 			if(prob->y[i] > 0)
// 			{
// 				if(fabs(alpha[i]) >= si.upper_bound_p)
// 					++nBSV;
// 			}
// 			else
// 			{
// 				if(fabs(alpha[i]) >= si.upper_bound_n)
// 					++nBSV;
// 			}
// 		}
// 	}

// 	info("nSV = %d, nBSV = %d\n",nSV,nBSV);

	decision_function f;
	f.alpha = alpha;
	f.rho = si.rho;
	return f;
}

//
// svm_model
//
struct svm_model
{
	svm_parameter param;	// parameter
	int nr_class;		// number of classes, = 2 in regression/one class svm
	int l;			// total #SV
	svm_node **SV;		// SVs (SV[l])
	double **sv_coef;	// coefficients for SVs in decision functions (sv_coef[n-1][l])
	double *rho;		// constants in decision functions (rho[n*(n-1)/2])

	// for classification only

	int *label;		// label of each class (label[n])
	int *nSV;		// number of SVs for each class (nSV[n])
				// nSV[0] + nSV[1] + ... + nSV[n-1] = l
	// XXX
	int free_sv;		// 1 if svm_model is created by svm_load_model
				// 0 if svm_model is created by svm_train
};



const char *svm_check_parameter(const svm_problem *prob, const svm_parameter *param)
{
	// svm_type

	int svm_type = param->svm_type;
	if(svm_type != C_SVC &&
	   svm_type != NU_SVC &&
	   svm_type != ONE_CLASS &&
	   svm_type != EPSILON_SVR &&
	   svm_type != NU_SVR)
		return "unknown svm type";
	
	// kernel_type
	
	int kernel_type = param->kernel_type;
	if(kernel_type != LINEAR &&
	   kernel_type != POLY &&
	   kernel_type != RBF &&
	   kernel_type != SIGMOID &&
	   kernel_type != R && 
	   kernel_type != LAPLACE&&
	   kernel_type != BESSEL&&
	   kernel_type != ANOVA)
		return "unknown kernel type";

	// cache_size,eps,C,nu,p,shrinking

	if(param->cache_size <= 0)
		return "cache_size <= 0";

	if(param->eps <= 0)
		return "eps <= 0";

	if(svm_type == C_SVC ||
	   svm_type == EPSILON_SVR ||
	   svm_type == NU_SVR)
		if(param->C <= 0)
			return "C <= 0";

	if(svm_type == NU_SVC ||
	   svm_type == ONE_CLASS ||
	   svm_type == NU_SVR)
		if(param->nu < 0 || param->nu > 1)
			return "nu < 0 or nu > 1";

	if(svm_type == EPSILON_SVR)
		if(param->p < 0)
			return "p < 0";

	if(param->shrinking != 0 &&
	   param->shrinking != 1)
		return "shrinking != 0 and shrinking != 1";


	// check whether nu-svc is feasible
	
	if(svm_type == NU_SVC)
	{
		int l = prob->l;
		int max_nr_class = 16;
		int nr_class = 0;
		int *label = Malloc(int,max_nr_class);
		int *count = Malloc(int,max_nr_class);

		int i;
		for(i=0;i<l;i++)
		{
			int this_label = (int)prob->y[i];
			int j;
			for(j=0;j<nr_class;j++)
				if(this_label == label[j])
				{
					++count[j];
					break;
				}
			if(j == nr_class)
			{
				if(nr_class == max_nr_class)
				{
					max_nr_class *= 2;
					label = (int *)realloc(label,max_nr_class*sizeof(int));
					count = (int *)realloc(count,max_nr_class*sizeof(int));
				}
				label[nr_class] = this_label;
				count[nr_class] = 1;
				++nr_class;
			}
		}
	
		for(i=0;i<nr_class;i++)
		{
			int n1 = count[i];
			for(int j=i+1;j<nr_class;j++)
			{
				int n2 = count[j];
				if(param->nu*(n1+n2)/2 > min(n1,n2))
				{
					free(label);
					free(count);
					return "specified nu is infeasible";
				}
			}
		}
	}

	return NULL; 
}


extern "C" {

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>  

  struct svm_node ** sparsify (double *x, int r, int c)
  {
    struct svm_node** sparse;
    int         i, ii, count;
    
    sparse = (struct svm_node **) malloc (r * sizeof(struct svm_node *));
    for (i = 0; i < r; i++) {
      /* determine nr. of non-zero elements */
      for (count = ii = 0; ii < c; ii++)
	if (x[i * c + ii] != 0) count++;
      
      /* allocate memory for column elements */
      sparse[i] = (struct svm_node *) malloc ((count + 1) * sizeof(struct svm_node));
      
      /* set column elements */
      for (count = ii = 0; ii < c; ii++)
	if (x[i * c + ii] != 0) {
	  sparse[i][count].index = ii;
	  sparse[i][count].value = x[i * c + ii];
	  count++;
	}
      
      /* set termination element */
      sparse[i][count].index = -1;
    }
    
    return sparse;
  }
  

struct svm_node ** transsparse (double *x, int r, int *rowindex, int *colindex)
{
    struct svm_node** sparse;
    int i, ii, count = 0, nnz = 0;

    sparse = (struct svm_node **) malloc (r * sizeof(struct svm_node*));
    for (i = 0; i < r; i++) {
        /* allocate memory for column elements */
        nnz = rowindex[i+1] - rowindex[i];
        sparse[i] = (struct svm_node *) malloc ((nnz + 1) * sizeof(struct svm_node));

        /* set column elements */
        for (ii = 0; ii < nnz; ii++) {
            sparse[i][ii].index = colindex[count];
            sparse[i][ii].value = x[count];
            count++;
        }

        /* set termination element */
        sparse[i][ii].index = -1;
    }

    return sparse;

}



  void solve_smo(const svm_problem *prob, const svm_parameter* param,
		 double *alpha, Solver::SolutionInfo* si, double C, double *linear_term)
  {
    int l = prob->l;
    int i;

    switch(param->svm_type)
      {
      case C_SVC:
	{ double Cp,Cn;
	  double *minus_ones = new double[l];
	  schar *y = new schar[l];
	  for(i=0;i<l;i++)
	    {
	      alpha[i] = 0;
	      minus_ones[i] = -1;
	      if(prob->y[i] > 0) y[i] = +1; else y[i]=-1;
	    }
	  if(param->nr_weight > 0)
	    {
	      Cp = C*param->weight[0];
	      Cn = C*param->weight[1];
	    }
	  else 
	    Cp = Cn = C;
	  Solver s; //have to weight cost parameter for multiclass. problems 
	  s.Solve(l, SVC_Q(*prob,*param,y), minus_ones, y,
		  alpha, Cp, Cn, param->eps, si, param->shrinking);
	  delete[] minus_ones;
	  delete[] y;
	}
	break;
      case NU_SVC:
	{
	  schar *y = new schar[l];
	  double nu = param->nu;
	  double sum_pos = nu*l/2;
	  double sum_neg = nu*l/2;
	  for(i=0;i<l;i++)
	    if(prob->y[i]>0) 
	      { 
		y[i] = +1;  
		alpha[i] = min(1.0,sum_pos);
		sum_pos -= alpha[i];
	      }
	    else {
	      y[i] = -1;
	      alpha[i] = min(1.0,sum_neg);
	      sum_neg -= alpha[i];
	  }
	  double *zeros = new double[l];
	  for(i=0;i<l;i++)
	    zeros[i] = 0;
	  Solver_NU s;
	  s.Solve(l, SVC_Q(*prob,*param,y), zeros, y,
		  alpha, 1.0, 1.0, param->eps, si,  param->shrinking);
	  double r = si->r;
	  //info("C = %f\n",1/r);
	  for(i=0;i<l;i++)
	    alpha[i] *= y[i]/r;
	  si->rho /= r;
	  si->obj /= (r*r);
	  si->upper_bound_p = 1/r;
	  si->upper_bound_n = 1/r;
	  delete[] y;
	  delete[] zeros;
	}
	break;
      case ONE_CLASS:
	{
	  double *zeros = new double[l];
	  schar *ones = new schar[l];
	  int n = (int)(param->nu*l);	// # of alpha's at upper bound
	  // set initial alpha probably usefull for smo 
	  for(i=0;i<n;i++)
	    alpha[i] = 1;
	  alpha[n] = param->nu * l - n;
	  for(i=n+1;i<l;i++)
	    alpha[i] = 0;
	  for(i=0;i<l;i++)
	    {
	      zeros[i] = 0;
	      ones[i] = 1;
	    }
	  
	  Solver s;
	  s.Solve(l, ONE_CLASS_Q(*prob,*param), zeros, ones,
		  alpha, 1.0, 1.0, param->eps, si, param->shrinking);

	  delete[] zeros;
	  delete[] ones;
	}
	  break;
      case EPSILON_SVR:
	{
	  double *alpha2 = new double[2*l];
	  double *linear_term = new double[2*l];
	  schar *y = new schar[2*l];
	  
	  for(i=0;i<l;i++)
	    {
	      alpha2[i] = 0;
	      linear_term[i] = param->p - prob->y[i];
	      y[i] = 1;
	      alpha2[i+l] = 0;
	      linear_term[i+l] = param->p + prob->y[i];
	      y[i+l] = -1;
	    }
	  Solver s;
	  s.Solve(2*l, SVR_Q(*prob,*param), linear_term, y,
		  alpha2, param->C, param->C, param->eps, si, param->shrinking);
	  double sum_alpha = 0;
	  for(i=0;i<l;i++)
	    {
	      alpha[i] = alpha2[i] - alpha2[i+l];
	      sum_alpha += fabs(alpha[i]);
	    }
	  //info("nu = %f\n",sum_alpha/(param->C*l));
	  
	  delete[] alpha2;
	  delete[] linear_term;
	  delete[] y; 
	}
	break;
		case NU_SVR:
		  {
		    double C = param->C;
		    double *alpha2 = new double[2*l];
		    double *linear_term = new double[2*l];
		    schar *y = new schar[2*l];
		    double sum = C * param->nu * l / 2;
		    for(i=0;i<l;i++)
		      {
			alpha2[i] = alpha2[i+l] = min(sum,C);
			sum -= alpha2[i];
			
			linear_term[i] = - prob->y[i];
			y[i] = 1;
			
			linear_term[i+l] = prob->y[i];
			y[i+l] = -1;
		      }
		    
		    Solver_NU s;
		    s.Solve(2*l, SVR_Q(*prob,*param), linear_term, y,
			    alpha2, C, C, param->eps, si, param->shrinking);
		    
		    //info("epsilon = %f\n",-si->r);
		    
		    for(i=0;i<l;i++)
		      alpha[i] = alpha2[i] - alpha2[i+l];
		    
		    delete[] alpha2;
		    delete[] linear_term;
		    delete[] y;
		  }
			break;
      }
  }

  SEXP smo_optim(SEXP x,
		 SEXP r, 
		 SEXP c, 
		 SEXP y, 
		 SEXP rowindex,
		 SEXP colindex,
		 SEXP sparse,
		 SEXP linear_term, 
		 SEXP kernel_type, 
		 SEXP svm_type, 
		 SEXP cost, 
		 SEXP nu, 
		 SEXP eps, 
		 SEXP gamma, 
		 SEXP degree, 
		 SEXP coef0, 
		 SEXP weightlabels, 
		 SEXP weights, 
		 SEXP nweights, 
		 SEXP cache,
		 SEXP epsilon, 
		 SEXP shrinking,
		 SEXP expr,
		 SEXP rho
		 )
  {
    
    SEXP alpha;
    double *alpha2;
    struct svm_parameter param;
    struct svm_problem  prob;
    int i;  
    const char* s;
    struct Solver::SolutionInfo si;
    param.svm_type    = *INTEGER(svm_type);
    param.kernel_type = *INTEGER(kernel_type); 
    param.degree      = *REAL(degree); 
    param.gamma       = *REAL(gamma);
    param.coef0       = *REAL(coef0);
    param.cache_size  = *REAL(cache);
    param.eps         = *REAL(epsilon);
    param.C           = *REAL(cost);
    param.nu          = *REAL(nu);
    param.expr        = expr;
    param.rho         = rho;
    param.nr_weight   = *INTEGER(nweights);
    if (param.nr_weight > 0) {
      param.weight      = (double *) malloc (sizeof(double) * param.nr_weight);
      memcpy (param.weight, REAL(weights), param.nr_weight * sizeof(double));
      param.weight_label = (int *) malloc (sizeof(int) * param.nr_weight);
      memcpy (param.weight_label, INTEGER(weightlabels), param.nr_weight * sizeof(int));
    }
    param.p           = *REAL(eps);
    param.shrinking   = *INTEGER(shrinking);
    param.lim = 1/(tgamma(param.degree+1)*pow(2,param.degree));    
    /* set problem */
    prob.l = *INTEGER(r);
    param.d = *INTEGER(c);
    // prob.y =  (double *) malloc (sizeof(double) * prob.l);
    prob.y = REAL(y);

    if (*INTEGER(sparse) > 0)
      prob.x = transsparse(REAL(x), *INTEGER(r), INTEGER(rowindex), INTEGER(colindex));
    else
      prob.x = sparsify(REAL(x), *INTEGER(r), *INTEGER(c)); 

    alpha2      = (double *) malloc (sizeof(double) * prob.l);
    s = svm_check_parameter(&prob, &param);
    if (s) {
      printf("%s",s);
    } else {

            solve_smo(&prob, &param, alpha2, &si, *REAL(cost), REAL(linear_term));
    }
    PROTECT(alpha = allocVector(REALSXP, prob.l+1));
    /* clean up memory */
    if (param.nr_weight > 0) {
      free(param.weight);
      free(param.weight_label);
    }
    for (i = 0; i < prob.l; i++) {free (prob.x[i]); REAL(alpha)[i] = *(alpha2+i); } 
    free (prob.x);
    REAL(alpha)[prob.l]=si.rho;
    free(alpha2); 
    UNPROTECT(1);  
    return alpha;
  }
}
