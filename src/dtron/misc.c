#include <math.h>
#include <stdlib.h>

void *xmalloc(size_t size)
{
        void *ptr = (void *) malloc(size);
        return ptr;
}
double mymax(double a, double b)
{
	if (a > b)
		return a;
	return b;
}
double mymin(double a, double b)
{
	if (a < b)
		return a;
	return b;	
}
double sign(double a, double b)
{
	if (b >= 0)
		return fabs(a);
	return -fabs(a);
}
