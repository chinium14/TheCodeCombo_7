/* ======================================================================== */
/* cwall.cpp                                                              */
/*            Driver for force and energy calculation with a smooth wall.   */
/* Adapted from code written by Tushar Jain.                                */
/* Added by Thomas Knotts 10 Mar 04.                                        */
/* ======================================================================== */
#ifdef WALL
#ifdef WWALL
#include "defines.h"

  /* ====================================================== */
  /* This calculates a purely repulsive wall potential.     */
  /* ====================================================== */


double cwall(int ibox){
  double ewall=0;
  int k = ibox;
  npr_1=0.81-0.008;//+0.01;
  npr_2=0.768;//-0.04; 
  npr_3=0.20 + 0.01;// - 0.1;// - 0.08;
  npr_w=0.01+0.005;//-0.002;// + 0.01;
  npr_p=0.05;
  double new_scale = 1.08;// 0.92;
  #if defined(WWALL)
  //printf("npr_1 = %f, npr_2 = %f, npr_3 = %f, npr_w = %f, npr_p = %f\n", npr_1, npr_2, npr_3, npr_w, npr_p); 
    double z_wall = wall[k].z[0];
    
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
        double sig = (sig_wall + pott[k][b].sig)/2.0;
//        printf("pott[k][b].sig = %f \n", pott[k][b].sig);
        double eps = sqrt(eps_wall * pott[k][b].eps);
//        printf("pott[k][b].eps = %f \n", pott[k][b].eps);
//        printf("eps = sqrt(eps_wall: %f * pott.eps: %f)",eps_wall, pott[k][b].eps);
//        double rc   = pow(2.0,(1.0/6.0))*sig;

//        double rc   = 25.0*sig;
//        double rc2  = rc * rc;
        double rho = wall[k].num_density[0][0];

//        double factor = 4.0*PI*eps*rho*pow(sig,3.0);
//        double factor2=factor/sig;
        double factor = PI*eps*rho*pow(sig,3.0)*new_scale;
        double factor2=factor/sig;
	
      #endif
      
      double dz = atom[k][b].z - z_wall;
//      printf("sig = %f, atom[k][b].z = %f, dz = %f \n",sig, atom[k][b].z, dz);
      double z_unit = dz / fabs(dz);
      double dr2 = dz*dz;
//      #ifdef WWALL
//      if(dr2 > rc2) continue;
//      #endif
      double sr2  = (sig * sig) / dr2;
      double sr3  = sr2 * sig / fabs(dz);
      double sr4  = sr2 * sr2;
      double sr6  = sr4 * sr2;
      double sr10 = sr6 * sr4;
      double sr9  = sr6 * sr3;
      double sr7  = sr4 * sr3;
      double sr8  = sr4 * sr4;
      double force;
//      printf("sr3 = %lf, sr7 = %lf, sr9 = %lf \n", sr3, sr7, sr9);
      #ifdef WWALL
     // printf("bgo[k].hydro_index[%d] is %lf \n", b, bgo[k].hydro_index[b]);
/*      double nu_nu = wall[k].hydro_index[0][0] * bgo[k].hydro_index[b];
      double nu_nu_sqrt;
      if (nu_nu >= 0){	
		nu_nu_sqrt = sqrt(nu_nu);
		}else{
		nu_nu_sqrt = -sqrt(-nu_nu);
		}
*/
/*      double shift = fabs(factor*(npr_r*1/45.0*pow((sig/rc),9.0) - (1.0/6.0 - npr_a * 5.0/12.0*nu_nu_sqrt)*pow((sig/rc),3.0)));
      ewall += factor*(npr_r*sr9/45.0 - (1.0/6.0 - npr_a * 5.0/12.0*nu_nu_sqrt)*sr3) + shift;
      force =  factor2*(npr_r*0.2*sr10 - (0.5 - npr_a * 5.0/4.0*nu_nu_sqrt)*sr4);
*/      
/*      double shift = fabs(factor*(npr_1*(13.0/45.0)*pow((sig/rc),9.0) - (npr_a*(9.0/14.0) + npr_h *nu_nu_sqrt)*pow((sig/rc),7.0) + npr_2*(2.0/3.0)*pow((sig/rc),3.0)));
      ewall += factor*(npr_1*(13.0/45.0)*sr9 - (npr_a*(9.0/14.0) + npr_h*nu_nu_sqrt)*sr7 + npr_2*(2.0/3.0)*sr3) + shift;
      force =  factor2*(npr_1*(13.0/5.0)*sr10 - (npr_a*(9.0/2.0) + npr_h*7.0*nu_nu_sqrt)*sr8 + 2.0*npr_2*sr4);
*/

/*Updated from above: 
	1) Change "+ shift" to "- shift", because the pmf curve from the new surface will be a little bit over 0 at a long distance. 
	2) Change nu_nu_sqrt to wall[k].hydro_index[0][0], in which way, the hydrophobicity of the surface is the only consideration.

	double shift = fabs(factor*(npr_1*(13.0/45.0)*pow((sig/rc),9.0) - (npr_a*(9.0/14.0) + npr_h * wall[k].hydro_index[0][0])*pow((sig/rc),7.0) + npr_2*(2.0/3.0)*pow((sig/rc),3.0)));
               ewall += factor*(npr_1*(13.0/45.0)*sr9 - (npr_a*(9.0/14.0) + npr_h * wall[k].hydro_index[0][0])*sr7 + npr_2*(2.0/3.0)*sr3) - shift;
               force =  factor2*(npr_1*(13.0/5.0)*sr10 - (npr_a*(9.0/2.0) + npr_h * 7.0 * wall[k].hydro_index[0][0])*sr8 + 2.0*npr_2*sr4);
*/

/*Updated from above:
	npr_a to npr_2, npr_2 to npr_3   add in npr_w, npr_p.
	wall[].hydro_index from 7th power term to 3rd power term.
	add in -4.5 for wall[].hydro[][] 
*/
//        printf("wall.hydro_index = %f, and bgo.hydro_index[b] = %f \n", wall[k].hydro_index[0][0], bgo[k].hydro_index[b]);
//	double third_term_co = npr_3*(2.0/3.0)-npr_w*(wall[k].hydro_index[0][0]-4.5)-npr_p*(bgo[k].hydro_index[b]);
//        printf("The coefficient of the  third power term = %f \n", third_term_co);
//        printf("The shift of the third power term = %f \n", pow((sig/rc),3.0));



//double shift = factor*(npr_1*(13.0/45.0)*pow((sig/rc),9.0) - npr_2*(9.0/14.0)*pow((sig/rc),7.0) + (npr_3*(2.0/3.0)-npr_w*(wall[k].hydro_index[0][0]-4.5)-npr_p*(bgo[k].hydro_index[b]))*pow((sig/rc),3.0));
       ewall+= factor*(npr_1*(13.0/45.0)*sr9 - npr_2*(9.0/14.0)*sr7 + (npr_3*(2.0/3.0)-npr_w*(wall[k].hydro_index[0][0]-4.5)-npr_p*(bgo[k].hydro_index[b]))*sr3);// - shift;
       force = factor2*(npr_1*(13.0/5.0)*sr10 - npr_2*(9.0/2.0)*sr8 + (npr_3*2.0-npr_w*3.0*(wall[k].hydro_index[0][0]-4.5)-npr_p*3.0*(bgo[k].hydro_index[b]))*sr4);

//      if (npr_2*(9.0/2.0)*sr8 < (npr_3*2.0-npr_w*3.0*(wall[k].hydro_index[0][0]-4.5)-npr_p*3.0*(bgo[k].hydro_index[b]))*sr4 ) 
/*
                  printf("b = %d, sig = %f, dz = %f, f4 = %f, f8 = %f, f10 = %f \n", 
                                                           b,
                                                           sig,
                                                           dz,
                                                           (npr_3*2.0-npr_w*3.0*(wall[k].hydro_index[0][0]-4.5)-npr_p*3.0*(bgo[k].hydro_index[b]))*sr4,
                                                           -npr_2*(9.0/2.0)*sr8, 
                                                           npr_1*(13.0/5.0)*sr10);
*/
/*Updated from above:
	move the npr_w part from the 3rd power term to 7th power term.
*/
/*	double shift = fabs(factor*(npr_1*(13.0/45.0)*pow((sig/rc),9.0) - (npr_2*(9.0/14.0) + npr_w * (wall[k].hydro_index[0][0] - 4.5))*pow((sig/rc),7.0) + npr_3*(2.0/3.0) *pow((sig/rc),3.0) - npr_p*(bgo[k].hydro_index[b])*pow((sig/rc),3.0)));
               ewall += factor*(npr_1*(13.0/45.0)*sr9 - (npr_2*(9.0/14.0) + npr_w * (wall[k].hydro_index[0][0] - 4.5))*sr7 + npr_3*(2.0/3.0) * sr3 - npr_p*(bgo[k].hydro_index[b]) * sr3) - shift;
               force =  factor2*(npr_1*(13.0/5.0)*sr10 - (npr_2*(9.0/2.0) + 7.0 * npr_w * (wall[k].hydro_index[0][0] - 4.5))*sr8 + npr_3*2.0 * sr4 - npr_p * 3.0 * (bgo[k].hydro_index[b]) * sr4);
*/

      //printf("nu_nu = %lf \n", nu_nu);
      #endif
      
      ff[k][b].z  +=  force * z_unit;
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
