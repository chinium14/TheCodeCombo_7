#ifdef NMA
#define NRANSI
#include "nrutil.h"
#include "nr.h"
#ifdef MKL
#include <mkl.h>
#else
#include <math.h>
#endif
#define TOL 1e-5

int ncom;
double *pcom,*xicom,(*nrfunc)(double []);

void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []))
{
	double brent(double ax, double bx, double cx,
		double (*f)(double), double tol, double *xmin);
	double f1dim(double x);
	void f2dim(double x);
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
		double *fc, double (*func)(double));
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}

	double temp=0.0;
	int index=0;
	nrfunc=func;
	for (j=1;j<=n;j++) 
	{
		if(temp<fabs(xicom[j]))
		{
			temp=fabs(xicom[j]);
			index=j;
		}
	}

	for (j=1;j<=n;j++) 
		xicom[j]/=(temp);

	ax=0.00;
	xx=0.02;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);

	for (j=1;j<=n;j++) {
		xi[j] *= xmin/(temp);
		p[j] += xi[j];
	}

	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
}
#undef TOL
#undef NRANSI
#endif

