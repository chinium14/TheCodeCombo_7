/* ======================================================================== */
/* vinit.cpp                                                                */
/*                                                                          */
/*		This subroutine generates initial velocites from a gaussian         */
/* distribution whose mean corresponds to the simulation temperature. It    */
/* is called if no simul.vel# is found in the INPUT folder. The net linear  */
/* momentum is also set to zero.											*/
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
double dgauss(double,double);
void zero_angular_momentum(int);
void pbc_all(int);

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void vinit (void) 
{


  /* ================================================================== */
  /*                                                                    */
  /* Initialize velocites from a Gaussing destribution.                 */
  /*                                                                    */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++) {
    for(int i=0; i<box[k].boxns; i++) {
      double imass = 1.0/pott[k][i].mas;
	  double ave = 0.0;
      double var = sqrt(RG*sim.T[k]*0.1*imass);
      vv[k][i].x = 1.0*dgauss(ave,var);
      vv[k][i].y = 1.0*dgauss(ave,var);
      vv[k][i].z = 1.0*dgauss(ave,var);
    }
  }


  for(int k=0; k<sim.NB; k++) {
  /* ----------------------------------------------- */
  /* If KONS is defined, set velocities of          */
  /* constrained atoms to zero.                      */
  /* ----------------------------------------------- */    
#ifdef KONS
	  for(int i=0; i<cons[k].n; i++){
		  int j=cons[k].site[i];
		  vv[k][j].x = 0.0;
		  vv[k][j].y = 0.0;
		  vv[k][j].z = 0.0;
	  }
#endif

  /* ================================================================== */
  /* Zero out the angular momentum.                                     */
  /* ================================================================== */
    //#ifndef WALL
    zero_angular_momentum(k);
    //#endif
  /* ================================================================== */
  /* Determine the net linear momentum.                                 */
  /* ================================================================== */

	  
  /* ----------------------------------------------- */
  /* Zero out the momentum accumulators.             */
  /* ----------------------------------------------- */    
    double pvx = 0.0;
    double pvy = 0.0;
    double pvz = 0.0;
    
  /* ----------------------------------------------- */
  /* Accumulate the momentum in each direaction.     */
  /* ----------------------------------------------- */    
    for(int i=0; i<box[k].boxns; i++) {
      pvx += vv[k][i].x * pott[k][i].mas;
      pvy += vv[k][i].y * pott[k][i].mas;
      pvz += vv[k][i].z * pott[k][i].mas;
    }

  /* ----------------------------------------------- */
  /* Determine the average net momentum in each      */
  /* direction.                                      */
  /* ----------------------------------------------- */    
#ifndef KONS
    pvx /= (1.0*box[k].boxns);
    pvy /= (1.0*box[k].boxns);
    pvz /= (1.0*box[k].boxns);
#endif
#ifdef KONS
    pvx /= (1.0*(box[k].boxns-cons[k].n));
    pvy /= (1.0*(box[k].boxns-cons[k].n));
    pvz /= (1.0*(box[k].boxns-cons[k].n));

#endif

  /* ================================================================== */
  /*                                                                    */
  /* Determine the current kinetic energy and zero out the net          */
  /* linear momentum.                                                   */
  /*                                                                    */
  /* ================================================================== */
    double kinet = 0.0;
    for(int i=0; i<box[k].boxns; i++) {
      double imass = 1.0 / pott[k][i].mas;
      vv[k][i].x -= pvx *imass;
      vv[k][i].y -= pvy *imass;
      vv[k][i].z -= pvz *imass;
	  #ifdef SMD
		  if(i== steermd[k].site1){
			  vv[k][i].x = 0.0;
			  vv[k][i].y = 0.0;
			  vv[k][i].z = 0.0;
			  uu[k][i].x = 0.0;
			  uu[k][i].y = 0.0;
			  uu[k][i].z = 0.0;
		  }
	  #endif
	  #ifdef NSMD
		  if(i== steermd[k].site1){
			  vv[k][i].x = 0.0;
			  vv[k][i].y = 0.0;
			  vv[k][i].z = 0.0;
			  uu[k][i].x = 0.0;
			  uu[k][i].y = 0.0;
			  uu[k][i].z = 0.0;
		  }
	  #endif
	  #ifdef KONS
	/* The variable i is changed to l here to avoid repitition	*/
	  for(int l=0; l<cons[k].n; l++){
		  int j=cons[k].site[l];
		  vv[k][j].x = 0.0;
		  vv[k][j].y = 0.0;
		  vv[k][j].z = 0.0;
	  }
      #endif

      kinet += pott[k][i].mas * (vv[k][i].x * vv[k][i].x +
				 vv[k][i].y * vv[k][i].y +
				 vv[k][i].z * vv[k][i].z);
    }

    kinet *= 0.5;

//    double nfree = (3.0 * box[k].boxns);
	double nfree = box[k].nfree;
#ifdef SMD
	nfree = nfree -3.0;
#endif
#ifdef NSMD
	nfree = nfree -6.0;
#endif

  /* ================================================================== */
  /*                                                                    */
  /* Determine the scale for the velocites to get the desired average   */
  /* temperature and scale the temperatures.                            */
  /*                                                                    */
  /* ================================================================== */
    double kin0 = sim.T[k] * nfree * 0.5 * RG * 0.10;
    double scal = sqrt(kin0 / kinet);

#ifndef SMD
#ifndef NSMD
    for(int i=0; i<box[k].boxns; i++) {
      vv[k][i].x *= scal;
      vv[k][i].y *= scal;
      vv[k][i].z *= scal;
    }
#endif//NO NSMD
#endif// NO SMD
#ifdef SMD
    for(int i=0; i<box[k].boxns; i++) {
      if(i!= steermd[k].site1){
		  vv[k][i].x *= scal;
		  vv[k][i].y *= scal;
		  vv[k][i].z *= scal;
	  }
    }
#endif//SMD
#ifdef NSMD
    for(int i=0; i<box[k].boxns; i++) {
      if(i!= steermd[k].site1 && i!= steermd[k].site2){
		  vv[k][i].x *= scal;
		  vv[k][i].y *= scal;
		  vv[k][i].z *= scal;
	  }
	  else if(i== steermd[k].site2){
			  vv[k][i].x = steermd[k].v*steermd[k].ex;
			  vv[k][i].y = steermd[k].v*steermd[k].ey;
			  vv[k][i].z = steermd[k].v*steermd[k].ez;
			  uu[k][i].x = vv[k][i].x;
			  uu[k][i].y = vv[k][i].y;
			  uu[k][i].z = vv[k][i].z;
	  }
    }
#endif//NSMD
  /* ================================================================== */
  /*                                                                    */
  /* Check the linear momentum again for debuggin purposes.             */
  /*                                                                    */
  /* ================================================================== */
/*    pvx = 0.0;
    pvy = 0.0;
    pvz = 0.0;
    
    for(int i=0; i<box[k].boxns; i++) {
      pvx += vv[k][i].x * pott[k][i].mas;
      pvy += vv[k][i].y * pott[k][i].mas;
      pvz += vv[k][i].z * pott[k][i].mas;
    }

    pvx /= (1.0*box[k].boxns);
    pvy /= (1.0*box[k].boxns);
    pvz /= (1.0*box[k].boxns);
    //double pva = sqrt(pvx*pvx + pvy*pvy + pvz*pvz);

    //    printf("%e  %e  %e  %e\n",pvx,pvy,pvz,pva);
*/
  /* ================================================================== */
  /*                                                                    */
  /* Set both velocity arrays to the new values.                        */
  /*                                                                    */
  /* ================================================================== */
    for(int i=0; i<box[k].boxns; i++) {
      uu[k][i].x = vv[k][i].x;
      uu[k][i].y = vv[k][i].y;
      uu[k][i].z = vv[k][i].z;
    }
#ifdef PR_NPT
	for(int i=0; i<box[k].boxns; i++){

	
		double mass = pott[k][i].mas;

		vvab[k][0] = vvab[k][0] + vv[k][i].x * vv[k][i].x * mass;
		vvab[k][1] = vvab[k][1] + vv[k][i].x * vv[k][i].y * mass;
		vvab[k][2] = vvab[k][2] + vv[k][i].x * vv[k][i].z * mass;
		vvab[k][4] = vvab[k][4] + vv[k][i].y * vv[k][i].y * mass;
		vvab[k][5] = vvab[k][5] + vv[k][i].y * vv[k][i].z * mass;
		vvab[k][8] = vvab[k][8] + vv[k][i].z * vv[k][i].z * mass;
	}
  /* ----------------------------------------------- */
  /* Convert gm/mol to kg.                           */
  /* ----------------------------------------------- */
    vvab[k][0] = vvab[k][0] / 1000.0 / NA;
		vvab[k][1] = vvab[k][1] / 1000.0 / NA;
		vvab[k][2] = vvab[k][2] / 1000.0 / NA;
		vvab[k][4] = vvab[k][4] / 1000.0 / NA;
		vvab[k][5] = vvab[k][5] / 1000.0 / NA;
		vvab[k][8] = vvab[k][8] / 1000.0 / NA;



	//to have a symmetric vvab matrix;
	vvab[k][3] = vvab[k][1];
	vvab[k][6] = vvab[k][2];
	vvab[k][7] = vvab[k][5];

#endif
  }

}

void zero_angular_momentum(int ibox){
  int k = ibox;
  double Lx, Ly, Lz;
  double A[3], B[3], C[3];
  double I[3][3] = {{0., 0., 0.},{0. ,0. ,0.},{0. ,0., 0.}}; 
  double II[3][3] = {{0., 0., 0.},{0. ,0. ,0.},{0. ,0., 0.}}; 
  double Omega[3] = {0.,0.,0.};
  int ns = box[k].boxns;

  /* ================================================================== */
  /* The linear momentum and COM must be zero before the angular        */
  /* momentum can be zeroed out.                                        */
  /* ================================================================== */
  /* ----------------------------------------------- */
  /* Zero out the momentum and COM accumulators.     */
  /* ----------------------------------------------- */    
    double pvx = 0.0;
    double pvy = 0.0;
    double pvz = 0.0;
 	  double xc = 0.0;
	  double yc = 0.0;
	  double zc = 0.0;
    double xco = 0.0;
    double yco = 0.0;
    double zco = 0.0;
	  double mc = 0.0;
    
  /* ----------------------------------------------- */
  /* Accumulate the momentum and COM in each         */
  /* direaction.                                     */
  /* ----------------------------------------------- */    
    for(int i=0; i<box[k].boxns; i++) {
      double mass = pott[k][i].mas;
      pvx += vv[k][i].x * mass;
      pvy += vv[k][i].y * mass;
      pvz += vv[k][i].z * mass;
		  xc += atnopbc[k][i].x * mass;
		  yc += atnopbc[k][i].y * mass;
		  zc += atnopbc[k][i].z * mass;
		  mc += pott[k][i].mas;
    }

  /* ----------------------------------------------- */
  /* Determine the average net momentum in each      */
  /* direction.                                      */
  /* ----------------------------------------------- */    
    pvx /= (1.0*box[k].boxns);
    pvy /= (1.0*box[k].boxns);
    pvz /= (1.0*box[k].boxns);
	  xc /= mc;
	  yc /= mc;
	  zc /= mc;
    xco = xc; //Store how much you moved in the x,y,z direction.
    yco = yc;
    zco = zc;

  /* ================================================================== */
  /* Zero out the net linear momentum and COM.                          */
  /* ================================================================== */
    for(int i=0; i<box[k].boxns; i++) {
      double imass = 1.0 / pott[k][i].mas;
      vv[k][i].x -= pvx *imass;
      vv[k][i].y -= pvy *imass;
      vv[k][i].z -= pvz *imass;
		  atnopbc[k][i].x -= xc;
		  atnopbc[k][i].y -= yc;
		  atnopbc[k][i].z -= zc;
		  atom[k][i].x     = atnopbc[k][i].x;
		  atom[k][i].y     = atnopbc[k][i].y;
		  atom[k][i].z     = atnopbc[k][i].z;
    }
	  pbc_all(k);
  /* ================================================================== */
  /* Check the linear momentum and COM again for debuggin purposes.     */
  /* ================================================================== */
/*    pvx = 0.0;
    pvy = 0.0;
    pvz = 0.0;
    
    for(int i=0; i<box[k].boxns; i++) {
      pvx += vv[k][i].x * pott[k][i].mas;
      pvy += vv[k][i].y * pott[k][i].mas;
      pvz += vv[k][i].z * pott[k][i].mas;
    }

    pvx /= (1.0*box[k].boxns);
    pvy /= (1.0*box[k].boxns);
    pvz /= (1.0*box[k].boxns);
    double pva = sqrt(pvx*pvx + pvy*pvy + pvz*pvz);


	  xc =0.0; yc =0.0; zc =0.0;
	  for(int i=0; i<box[k].boxns; i++) {
		  xc += atnopbc[k][i].x * pott[k][i].mas;
		  yc += atnopbc[k][i].y * pott[k][i].mas;
		  zc += atnopbc[k][i].z * pott[k][i].mas;
		  mc += pott[k][i].mas;
	  }
*/


  /* ================================================================== */
  /* Now correct the angular momentum.                                  */
  /* ================================================================== */
    //for(int j = 0; j<2; j++){
  Lx = 0.;
 	Ly = 0.;
	Lz = 0.;
 	for(int i=0; i<ns; i++)
 	{
    double mass = pott[k][i].mas;
    double rx   = atom[k][i].x;
    double ry   = atom[k][i].y;
    double rz   = atom[k][i].z;
  	//Calculate the total angular momentum
    A[0] = rx;
		A[1] = ry;
		A[2] = rz;
		B[0] = mass * vv[k][i].x; 
		B[1] = mass * vv[k][i].y;  
		B[2] = mass * vv[k][i].z;  
		CROSS(A,B,C);
		Lx += C[0];
		Ly += C[1];
		Lz += C[2];
    //Calculate the Inertia tensor
		I[0][0] += mass * (ry * ry + rz * rz); 
		I[1][1] += mass * (rx * rx + rz * rz); 
		I[2][2] += mass * (rx * rx + ry * ry); 
		I[0][1] -= mass * (rx * ry); 
		I[0][2] -= mass * (rx * rz); 
    I[1][2] -= mass * (ry * rz); 
  }
  I[1][0] = I[0][1];
  I[2][0] = I[0][2];
  I[2][1] = I[1][2];
    
 	//Take inverse of the matrix 
 	INV(I,II);

 	//Calculate Omega
 	Omega[0] = II[0][0] * Lx + II[0][1] * Ly + II[0][2] * Lz;
 	Omega[1] = II[1][0] * Lx + II[1][1] * Ly + II[1][2] * Lz;
 	Omega[2] = II[2][0] * Lx + II[2][1] * Ly + II[2][2] * Lz;
 	for(int i=0; i<ns; i++)
 	{
		B[0] = atom[k][i].x;
		B[1] = atom[k][i].y;
		B[2] = atom[k][i].z;
		CROSS(Omega,B,C);
		vv[k][i].x -= C[0];
		vv[k][i].y -= C[1];
		vv[k][i].z -= C[2];
  }

  /* ================================================================== */
  /* Move the COM back to starting position.                            */
  /* ================================================================== */
  for(int i=0; i<ns; i++){
      atnopbc[k][i].x += xc;
      atnopbc[k][i].y += yc;
      atnopbc[k][i].z += zc;
      atom[k][i].x     = atnopbc[k][i].x;
      atom[k][i].y     = atnopbc[k][i].y;
      atom[k][i].z     = atnopbc[k][i].z;
  }
  pbc_all(k);


    //}
  /* Check to make sure angular momentum is zero */
/*  Lx = 0.;
  Ly = 0.;
  Lz = 0.;
 	for(int i=0; i<ns; i++)
 	{
    double mass = pott[k][i].mas;
    double rx   = atom[k][i].x;
    double ry   = atom[k][i].y;
    double rz   = atom[k][i].z;
  	//Calculate the total angular momentum
    A[0] = rx;
		A[1] = ry;
		A[2] = rz;
		B[0] = mass * vv[k][i].x; 
		B[1] = mass * vv[k][i].y;  
		B[2] = mass * vv[k][i].z;  
		CROSS(A,B,C);
		Lx += C[0];
		Ly += C[1];
		Lz += C[2];
  }
*/
}
