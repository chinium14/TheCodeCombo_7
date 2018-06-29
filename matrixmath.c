/* ======================================================================== */
/* matrixmath.cpp                                        					*/
/* This file contains subroutines to handle mathmatics of matrix.  jacobi() */
/* computes all the eigenvalues and vectors of a symmetric matrix.          */
/* ======================================================================== */

/* ======================================================================== */
/* eigen()                                                                  */
/* Computes all eigenvalues and eigenvectors of a real symmetric matrix     */
/* a[1..n][1..n]. On output, elements of a above the diagonal are           */
/* destroyed. d[1..n] returns the eigenvalues of a.  v[1..n][1..n] is a     */
/* matrix whose columns contain, on output, the normalized eigenvectors of  */
/* a. nrot returns the number of Jacobi rotations that were required.       */
/* Based on Numerical Recipies in C "jacobi" p467.                          */
/* ======================================================================== */
#ifdef PR_NPT
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);
#include <math.h>
void eigen(double a[][3],int n,double d[],double v[][3]){

	int j,iq,ip,i,nrot;
	double tresh,theta,tau,t,sm,s,h,g,c,b[3],z[3];
	
	n=n-1;

	for (ip=0;ip<=n;ip++){              /*           Initialize to the identity matrix.*/
		for (iq=0;iq<=n;iq++)
			v[ip][iq]=0.0;
			v[ip][ip]=1.0;
		}

	for (ip=0;ip<=n;ip++){              /* Initialize b and d to the diagonal of a*/
             b[ip]=d[ip]=a[ip][ip];
             z[ip]=0.0;
	}	
	nrot=0;

	for (i=1;i<=50;i++){
		sm=0.0;
		for (ip=0;ip<=n-1;ip++){                                                                                                                                  /* Sum of-diagonal elements.*/
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
			}
			if (sm == 0.0){
				b[0]=0;
				z[0]=0;
				b[1]=0;
				z[1]=0;
				b[2]=0;
				z[2]=0;                 /* The normal return, which relies on quadratic convergence to machine underflow.*/       
				return;
			}

			if (i < 4)
				tresh=0.2*sm/(n*n);
			else
				tresh=0.0;
			for (ip=0;ip<=n-1;ip++){
				for (iq=ip+1;iq<=n;iq++){
					g=100.0*fabs(a[ip][iq]);
					if(i > 4 &&(float)(fabs(d[ip])+g) == (float)fabs(d[ip])
					&& (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
                        a[ip][iq]=0.0;
					else
                           if (fabs(a[ip][iq]) > tresh) 
                           {
                            h=d[iq]-d[ip];
                            if ((float)(fabs(h)+g) == (float)fabs(h))
                                        t=(a[ip][iq])/h;
                            else
                            {
                                        theta=0.5*h/(a[ip][iq]);
                                        t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                                        if (theta < 0.0) t = -t;
                            }
							c=1.0/sqrt(1+t*t);
							s=t*c;
							tau=s/(1.0+c);
							h=t*a[ip][iq];
							z[ip] -= h;
							z[iq] += h;
							d[ip] -= h;
							d[iq] += h;
							a[ip][iq]=0.0;
							for (j=0;j<=ip-1;j++){ 
								ROTATE(a,j,ip,j,iq)
							}
							for (j=ip+1;j<=iq-1;j++){ 
								ROTATE(a,ip,j,j,iq)
							}
							for (j=iq+1;j<=n;j++){
								ROTATE(a,ip,j,iq,j)
							}
							for (j=0;j<=n;j++){
								ROTATE(v,j,ip,j,iq)
							}
							nrot++;
						   }
				}
			}
 
			for (ip=0;ip<=n;ip++){
				b[ip] += z[ip];
				d[ip]=b[ip];
				z[ip]=0.0;
			}
	}
}


/* ======================================================================== */
/* matrixmul()                                                              */
/* This subroutine multiplies two 3x3 matrices together.                    */
/* INPUT	: a1, a1 (matrices), l(option)                                  */
/* OUTPUT	: a3 = a1 * a2 (l=1), a3 = a1(transpose) * a2 (l=2)             */               
/* ======================================================================== */
void matrixmul(double a1[][3], double a2[][3], double a3[][3], int l){
	
	int i, j, k;
	double tmp = 0.0;

	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			for(k=0; k<3; k++){
				if(l==1) tmp = tmp + a1[i][k] * a2[k][j];
				else tmp = tmp + a1[k][i] * a2[k][j];
			}
			a3[i][j] = tmp;
			tmp = 0.0;
		}
	}
}
#endif

