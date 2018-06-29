#ifdef CONFIGT
/* ======================================================================== */
/* cwall.cpp                                                              */
/*            Driver for force and energy calculation with a smooth wall.   */
/* Adapted from code written by Tushar Jain.                                */
/* Added by Thomas Knotts 10 Mar 04.                                        */
/* ======================================================================== */
#ifdef WALL
#include "defines.h"

  /* ====================================================== */
  /* This calculates a purely repulsive wall potential.     */
  /* ====================================================== */


double cwall2(int ibox){
  double ewall=0;
  int k = ibox;
  #if defined(AWALL) || defined(NWALL) || defined(RWALL)
    double sig = 0.5*(golik[k].d_repulsive + wall[k].sig[0][0]);
    double eps = sqrt(golik[k].eps * wall[k].eps[0][0]);
    double rc   = pow((2.0/5.0),(1.0/6.0))*sig;
    double rc2  = rc * rc;
    double rho = wall[k].num_density[0][0];
    #ifdef RWALL
      double factor = 4.0*PI*eps*rho*pow(sig,3.0)/45.0;
      double factor2=factor/sig*9.0;
      double factor3=factor2/sig*10.0;
    #else
      double factor = 4.0*PI*eps*rho*pow(sig,3.0);
      double factor2=factor/sig;
      double factor3=factor2/sig*2.0;
    #endif
    #ifdef NWALL
      double shift = fabs(factor*(1/45.0*pow((sig/rc),9.0) - (1.0/6.0)*pow((sig/rc),3.0)));
    #endif
    double z_wall = wall[k].z[0];

    for(int b=0; b<box[k].boxns; b++){
	    double dz = atom[k][b].z - z_wall;
	    double z_unit = dz / fabs(dz);
      double dr2 = dz*dz;
      double dr = sqrt(dr2);
      #ifdef NWALL
        if(dr2 > rc2) continue;
      #endif
      double sr   = sig/dr;
		  double sr2  = (sig * sig) / dr2;
		  double sr3  = sr2 * sr;
		  double sr4  = sr2 * sr2;
		  double sr6  = sr4 * sr2;
		  double sr10 = sr6 * sr4;
		  double sr9  = sr6 * sr3;
      double sr11 = sr10* sr;
      double sr5  = sr4 * sr;
		  double force,hesz;
	  
		  #ifdef NWALL
		    ewall += factor*(sr9/45.0 - (1.0/6.0)*sr3) + shift;
 		    force =  factor2*(0.2*sr10 - 0.5*sr4);
        hesz = factor3*(sr11-sr5);
		  #elif defined(AWALL)
		    ewall += factor*(sr9/45.0 - (1.0/6.0)*sr3);
 		    force =  factor2*(0.2*sr10 - 0.5*sr4);
        hesz = factor3*(sr11-sr5);
		  #elif defined(RWALL)
		    ewall += factor*sr9;
		    force = factor2*sr10;
        hesz = factor3*sr11;
		  #endif

		  ff[k][b].z  +=  force * z_unit;
      config[k].hesz	+= hesz;

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
    return(ewall + ebend);

  #elif defined(HWALL)
    double z_wall = wall[k].z[0];
    for(int i=0; i<box[k].boxns; i++){
  	  if(atom[k][i].z < z_wall){
        return(ELIMIT);
      }
    }
    return(0.0);

  #else //CWALL
    double sig = wall[k].sig[0][0];
    double eps = wall[k].eps[0][0];
    double eps48 = 48.0*eps;
    double z_wall = wall[k].z[0];
    double Rs = fabs(z_wall);
    double rc = sig*pow(2.0,(1.0/6.0))+Rs;
    double rc2= rc*rc;
    double dra2i, dra, drai, dr, dri;
    double sr, sr2, sr4, sr6, sr12, force;
    double term1,term2,eps48dri;
    /* ====================================================== */
    /* Calculate the distance from each particle to the       */
    /* origin and apply the potential from a curved surface.  */
    /* ====================================================== */
    for(int b=0; b<box[k].boxns; b++){
      double dx   = atom[k][b].x;
      double dy   = atom[k][b].y;
      double dz   = atom[k][b].z-z_wall;
      double dx2  = dx*dx;
      double dy2  = dy*dy;
      double dz2  = dz*dz;
      double dra2 = dx2 + dy2 + dz2;

      int flag = int(dra2/rc2); //flag is 0 within cut-off

      switch(flag){
      case 0:
        dra2i= 1.0/dra2;
        dra  = sqrt(dra2);
        drai = 1.0/dra;
        dr   = dra - Rs;
        dri  = 1.0/dr;

        sr   = sig*dri;
        sr2  = sr*sr;
        sr4  = sr2*sr2;
        sr6  = sr4*sr2;
        sr12 = sr6*sr6;

        eps48dri = eps48*dri;
      

        force  = eps48dri*(sr12-0.5*sr6)*drai;
        ewall += 4.0*eps*(sr12-sr6) + eps;

        term1 = force + eps48dri*dri*(13.0*sr12-3.5*sr6);
        term2 = term1*dra2i;

        config[k].hesx += -force + term2*dx2;
        config[k].hesy += -force + term2*dy2;
        config[k].hesz += -force + term2*dz2;

        ff[k][b].x += force * dx;
        ff[k][b].y += force * dy;
        ff[k][b].z += force * dz;

        break;
      default:
        break;
      }
    
    }


    return(ewall);
  #endif
}

#endif
#endif
