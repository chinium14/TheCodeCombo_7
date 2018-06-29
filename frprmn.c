#ifdef NMA
#define NRANSI
#include "nrutil.h"
#define ITMAX 200000
#define EPS 1.0e-10
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);

#include <stdio.h>
#ifdef MKL
#include <mkl.h>
#else
#include <math.h>
#endif
void frprmn(double p[], int n, double ftol, int *iter, double *fret,double (*func)(double []), void (*dfunc)(double [], double []),int itera)
{
	void linmin(double p[], double xi[], int n, double *fret,double (*func)(double []));
	int j,its,cnt=0;
	double gg,gam,fp,dgg;
	double *g,*h,*xi;
	
	g=vector(1,n);
	h=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,xi);

	for (j=1;j<=n;j++) 
	{
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}

	for (its=1;its<=itera;its++)
	{
//		dsy_c(its);
		*iter=its;
		linmin(p,xi,n,fret,func);

		if(its%10==0)
			printf("%lf %lf\n",fp,*fret);

		if(fabs(*fret-fp)<1e-09 )
			cnt++;
		else
			cnt=0;

		if(cnt==10)
		{
			FREEALL
			return;
		}

		fp=*fret;
		(*dfunc)(p,xi);
		
		dgg=gg=0.0;

		for (j=1;j<=n;j++) 
		{
			gg += g[j]*g[j];
//			dgg+=xi[j]*xi[j];
			dgg += (xi[j]+g[j])*xi[j];
		}

		if (gg == 0.0) 
		{
			FREEALL
			printf("Exit 2");
			return;
		}

		gam=dgg/gg;
		for (j=1;j<=n;j++) 
		{
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	printf("Exit 3");
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI
#endif

