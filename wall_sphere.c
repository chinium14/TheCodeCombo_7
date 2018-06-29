/* ======================================================================== */
/* cwall.cpp                                                              */
/*            Driver for force and energy calculation with a smooth wall.   */
/* Adapted from code written by Tushar Jain.                                */
/* Added by Thomas Knotts 10 Mar 04.                                        */
/* ======================================================================== */
#ifdef WALL
#ifdef SPHERE
#include "defines.h"

  /* ====================================================== */
  /* This calculates a purely repulsive wall potential.     */
  /* ====================================================== */


double cwall(int ibox){
  double ewall=0;
  int k = ibox;
  npr_1=0.81-0.008;
  npr_2=0.768;
  npr_3=0.20+0.01;
  npr_w=0.01+0.005;
  npr_p=0.05;
  
  double new_scale = 1.08;
  #if defined(SPHERE)
    
//    double z_wall = wall[k].z[0];

  struct sphere s0 = {0.0, 0.0, -80.0, 80.0};
    #ifndef GOLIK //if GOLIK not defined, do mixing rules on the fly
    double sig_wall = wall[k].sig[0][0];
    double eps_wall = wall[k].eps[0][0];
    #endif
    #ifndef XWHAM
    for(int b=0; b<box[k].boxns; b++){
    #endif
    #ifdef XWHAM
    for(int b=0; b<box[k].boxns-1; b++){
    #endif
      #ifndef GOLIK
        double sig = (sig_wall + pott[k][b].sig)*0.5;
        double eps = sqrt(eps_wall * pott[k][b].eps);
//        double rc   = pow(2.0,(1.0/6.0))*sig;
//        double rc   = 25.0*sig;
//        double rc   = 25.0*sig;
//        double rc2  = rc * rc;
        double rho = wall[k].num_density[0][0];

//        double factor = 4.0*PI*eps*rho*pow(sig,3.0);
//        double factor2=factor/sig;
        double factor = PI*eps*rho*pow(sig,3.0)*new_scale;
        double factor2=factor/sig;
//        double factor = 2*PI*eps*rho*s0.r*pow(sig,6.0);
//        double factor2=factor/sig;
	
      #endif
      
        double xx = atom[k][b].x - s0.x;
        double yy = atom[k][b].y - s0.y;
        double zz = atom[k][b].z - s0.z;
/*
      double x_unit = xx / fabs(xx);
      double y_unit = yy / fabs(yy);
      double z_unit = zz / fabs(zz);
*/
        double r2 = xx*xx + yy*yy + zz*zz;
        double r = sqrt(r2);
        double r_sphere = s0.r; 
        double r_diff = fabs(r - r_sphere);
//      double rs2 = r_sphere * r_sphere;
        #ifdef SPHERE
//        if(r_diff > rc) continue;
        #endif
        double sr = sig / r_diff;
        double sr2  = sr * sr;
        double sr3  = sr2 * sr;
        double sr4  = sr2 * sr2;
//      double sr5  = sr2 * sr3; 
        double sr6  = sr4 * sr2;
        double sr7  = sr4 * sr3;
        double sr10 = sr6 * sr4;
        double sr9  = sr6 * sr3;
        double sr8  = sr4 * sr4;
/* 
      double srs  = 1.0 / (r2 - rs2);
      double srs2 = srs * srs;
      double srs3 = srs2 * srs;
      double srs4 = srs2 * srs2;
      double srs5 = srs4 * srs; 
      double srs6 = srs4 * srs2;
      double sig2 = sig * sig;
      double sig4 = sig2 * sig2;
      double sig6 = sig2 * sig4;
*/
        double force;

        #ifdef SPHERE
/*	double shift = fabs(factor*1.0/r*(npr_1*sig6*(13.0/10.0)*(pow((1.0/(rc-r_sphere)),10.0)-pow((1.0/(rc*rc - r_sphere*r_sphere)),5.0))
			- (npr_2*sig4*(9.0/4.0))*(pow((1.0/(rc-r_sphere)),8.0)-pow((1.0/(rc*rc - r_sphere*r_sphere)),4.0))
			+ (npr_3)*(pow((1.0/(rc-r_sphere)),4.0)-pow((1.0/(rc*rc - r_sphere*r_sphere)),2.0)) 
			- npr_w * (wall[k].hydro_index[0][0] - 4.5) *(pow((1.0/(rc-r_sphere)),4.0)-pow((1.0/(rc*rc - r_sphere*r_sphere)),2.0)) 
			- npr_p*(bgo[k].hydro_index[b])*(pow((1.0/(rc-r_sphere)),4.0)-pow((1.0/(rc*rc - r_sphere*r_sphere)),2.0))));


               ewall += factor*1.0/r*(npr_1*sig6*(13.0/10.0)*(sr10-srs5) 
				- (npr_2*sig4*(9.0/4.0))*(sr8-srs4)
				+ (npr_3) * (sr4-srs2)
				- npr_w * (wall[k].hydro_index[0][0] - 4.5) * (sr4-srs2) 
				- npr_p*(bgo[k].hydro_index[b]) * (sr4-srs2)) - shift;


               force =  factor*(npr_1*(13.0/10.0)*(-1.0/r2*(sr10-srs5)+10.0*(srs6 - 1.0/ r *sr11))
				- (npr_2*(9.0/4.0))*(-1.0/r2*(sr8-srs4)+8.0*(srs5 - 1.0/ r *sr9))
				+ (npr_3) * (-1.0/r2*(sr4-srs2)+4.0*(srs3 - 1.0/ r *sr5)) 
				- npr_w * (wall[k].hydro_index[0][0] - 4.5) * (-1.0/r2*(sr4-srs2)+4.0*(srs3 - 1.0/ r *sr5))  
				- npr_p * (bgo[k].hydro_index[b]) * (-1.0/r2*(sr4-srs2)+4.0*(srs3 - 1.0/ r *sr5)));
*/
/*Updated from above:
	move the npr_w part from the 3rd power term to 7th power term.
*/
		
//        double shift = fabs(factor*(r_sphere/(rc+r_sphere))*(npr_1*(13.0/45.0)*pow((sig/rc),9.0) - (npr_2*(9.0/14.0))*pow((sig/rc),7.0) + (npr_3*(2.0/3.0) - npr_w * (wall[k].hydro_index[0][0] - 4.5)) *pow((sig/rc),3.0) - npr_p*(bgo[k].hydro_index[b])*pow((sig/rc),3.0)));
               ewall += factor*(r_sphere/r)*(npr_1*(13.0/45.0)*sr9 - (npr_2*(9.0/14.0))*sr7 + (npr_3*(2.0/3.0) - npr_w * (wall[k].hydro_index[0][0] - 4.5)) * sr3 - npr_p*(bgo[k].hydro_index[b]) * sr3);// - shift;
               force =  factor2*(r_sphere/r)*(npr_1*(13.0/5.0)*sr10 - (npr_2*(9.0/2.0))*sr8 + (npr_3*2.0 - npr_w * 3.0 * (wall[k].hydro_index[0][0] - 4.5)) * sr4 - npr_p * 3.0 * (bgo[k].hydro_index[b]) * sr4);

      //printf("nu_nu = %lf \n", nu_nu);
        #endif
      
//      ff[k][b].x  +=  force * (xx / r) * x_unit;
//      ff[k][b].y  +=  force * (yy / r) * y_unit;
//      ff[k][b].z  +=  force * (zz / r) * z_unit;
        ff[k][b].x  +=  force * (xx / r);
        ff[k][b].y  +=  force * (yy / r);
        ff[k][b].z  +=  force * (zz / r);
    }

    double ebend = 0.0;
    if(wall[k].angle_k>0.0){
      /* ====================================================== */
      /* This calculates an angle potential to keep the         */
      /* the free site at a fixed angle from the surface.       */
      /* ====================================================== */
      int site = wall[k].angle_site;

      /* ----------------------------------------- */
      /* Calculate the distances between the atoms */
      /* involved in the bending.                  */
      /* ----------------------------------------- */
      double dx = atom[k][site].x;
      double dy = atom[k][site].y;
      double dz = atom[k][site].z;

      /* ----------------------------------------- */
      /* Calculate the angle of the bend and then  */
      /* the energy using the CHARMM force field.  */
      /* ----------------------------------------- */
      double x2y2 = dx * dx + dy * dy;
      double sqrt_x2y2 = sqrt(x2y2);
      double tanphi = dz/sqrt_x2y2;
      double phi = atan(tanphi);
      double cosphi = cos(phi);
      double cosphi2 = cosphi * cosphi;

      double rr = phi - wall[k].angle;
      double fr = 2.0*wall[k].angle_k * rr * cosphi2 * dz / sqrt_x2y2;
      ebend = wall[k].angle_k * rr * rr;

      /* ----------------------------------------- */
      /* Calculate the forces.                     */
      /* ----------------------------------------- */
      ff[k][site].x += fr * dx / x2y2;
      ff[k][site].y += fr * dy / x2y2;
      ff[k][site].z += -fr / dz;
    }
//    return(ewall + ebend);
    return(ewall);
  #endif
}
#endif
#endif
