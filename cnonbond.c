/* ======================================================================== */
/* cnonbond.cpp                                                             */
/*                                                                          */
/*		This subroutine calculates the nonbonded stretching energies and    */
/* forces. This includes the van der Waals attraction and repulsion and the */
/* colombic interactions.  It returns the no value, but stores the energies */
/* in an array "energy" which is defined in forces() or force_long().  It   */
/* is called from either forces() or force_long(). This subroutine is       */
/* called if NLIST is not defined.  If NLIST is defined, cnonbond_nblist()  */
/* is called.*/
/*                                                                          */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                      *energy:	A pointer to an array defined in        */
/*                                  forces() or force_long().  It contains  */
/*                                  three entries.  Two for unshifted and   */
/*                                  shifted vdw energies and one for        */
/*                                  coulomib energies.                      */
/*                                                                          */
/* ======================================================================== */

#ifndef NLIST

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void cnonbond (int ibox, double *energy)
{

  int k = ibox;

  /* ================================================================== */
  /*                                                                    */
  /* Zero out the energy accumulator                                    */
  /*                                                                    */
  /* ================================================================== */
  double enonbond  = 0.0;
  double enonbonds = 0.0;
  double ecoulomb  = 0.0;

  /* ================================================================== */
  /*                                                                    */
  /* Zero out the virial accumulators                                   */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  double wxx = 0.0;
  double wyx = 0.0;
  double wyy = 0.0;
  double wzx = 0.0;
  double wzy = 0.0;
  double wzz = 0.0;
#endif

  /* ================================================================== */
  /*                                                                    */
  /* Get and assign information for pbc and cutoffs                     */
  /*                                                                    */
  /* ================================================================== */
  double rc   = sqrt(box[k].rc2);
  double rc2  = box[k].rc2;
  double rfcs = (sim.epsRF[k] - 1.0) / (2.0*sim.epsRF[k] + 1.0);
  double rfc  = rfcs / (rc2*rc);
  double bl   = box[k].boxl;
  double bh   = box[k].boxh;
 
  /* ================================================================== */
  /*                                                                    */
  /* Calculate the energy and force associated with each pair.  The     */
  /* loops cover each pair only once.                                   */
  /*                                                                    */
  /* ================================================================== */
  for(int a=0; a<box[k].boxns-1; a++) {

  /* ----------------------------------------------- */
  /* Zero out the force accumulators for particle a. */
  /* ----------------------------------------------- */
	double fxa = 0.0;
    double fya = 0.0;
    double fza = 0.0;
  
  /* ----------------------------------------------- */
  /* Assign position of particle a.                  */
  /* ----------------------------------------------- */    
    double ax = atom[k][a].x;
    double ay = atom[k][a].y;
    double az = atom[k][a].z;

  /* ----------------------------------------------- */
  /* Loop around atoms "above" particle a in array.  */
  /* ----------------------------------------------- */
    for(int b=a+1; b<box[k].boxns; b++) {

      double fx, fy, fz;

  /* ----------------------------------------------- */
  /* Calculate the x,y,z distances between           */
  /* particles a & b.                                */
  /* ----------------------------------------------- */
      double dx = ax - atom[k][b].x;
      double dy = ay - atom[k][b].y;
      double dz = az - atom[k][b].z;

  /* ----------------------------------------------- */
  /* Apply minimum image convention.                 */
  /* ----------------------------------------------- */
      double pbcx = 0.0;
      double pbcy = 0.0;
      double pbcz = 0.0;

      if(dx >  bh) pbcx =- bl;
      if(dx < -bh) pbcx =+ bl;
      if(dy >  bh) pbcy =- bl;
      if(dy < -bh) pbcy =+ bl;
      if(dz >  bh) pbcz =- bl;
      if(dz < -bh) pbcz =+ bl;

      dx += pbcx;
      dy += pbcy;
      dz += pbcz;

  /* ----------------------------------------------- */
  /* Calculate the distance between particles a & b. */
  /* If the distance is less than the cutoff,        */
  /* calculated the energy and the force. If it is   */
  /* greater, go to the next pair in the loop.       */
  /* ----------------------------------------------- */
      double dr2 = dx*dx + dy*dy + dz*dz;
      
      if(dr2 > rc2) continue;

      printf("%d %d  %lf \n",a,b,sqrt(dr2));
      
      double eps = ljset[k].pot[a][b].eps;
      double sig = ljset[k].pot[a][b].sig;
      double qq  = ljset[k].pot[a][b].qq;		  

      double sr2  = (sig * sig) / dr2;
      double sr6  = sr2 * sr2 * sr2;
      double sr12 = sr6 * sr6;

      double sr2s  = (sig * sig) / rc2;
      double sr6s  = sr2s * sr2s * sr2s;
      double sr12s = sr6s * sr6s;

      double force = (12.0*eps/dr2) * (sr12 - sr6);
	  enonbond  +=  eps * (sr12 - 2*sr6);
      enonbonds +=  eps * (sr12s - 2*sr6s);
  /* ----------------------------------------------- */
  /* Calculate Coulombic interactions.               */
  /* ----------------------------------------------- */
#ifdef COULOMB
      double dri = 1.0 / sqrt(dr2);
      force += (qq/dr2) * (dri - 2.0 * rfc * dr2);
      enonbond  += qq * (dri + rfc * dr2);
	  ecoulomb += qq * (dri + rfc * dr2)  - qq * (1.0 + rfc) / rc;;
      enonbonds += qq * (1.0 + rfc) / rc; 
#endif


  /* ----------------------------------------- */
  /* Accumulate the forces.                    */
  /* ----------------------------------------- */
      fx          = force * dx;
      fxa        += fx;
      ff[k][b].x -= fx;

      fy          = force * dy;
      fya        += fy;
      ff[k][b].y -= fy;

      fz          = force * dz;
      fza        += fz;
      ff[k][b].z -= fz;

  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
      double sxx = fx * dx;
      double syx = fy * dx;
      double syy = fy * dy;
      double szx = fz * dx;
      double szy = fz * dy;
      double szz = fz * dz;
	wxx += sxx;
	wyx += syx;
	wyy += syy;
	wzx += szx;
	wzy += szy;
	wzz += szz;
#ifdef ZHOU
	stresses[k][a][0][0] += sxx * 0.5;
	stresses[k][a][1][0] += syx * 0.5;
	stresses[k][a][1][1] += syy * 0.5;
	stresses[k][a][2][0] += szx * 0.5;
	stresses[k][a][2][1] += szy * 0.5;
	stresses[k][a][2][2] += szz * 0.5;

	stresses[k][b][0][0] += sxx * 0.5;
	stresses[k][b][1][0] += syx * 0.5;
	stresses[k][b][1][1] += syy * 0.5;
	stresses[k][b][2][0] += szx * 0.5;
	stresses[k][b][2][1] += szy * 0.5;
	stresses[k][b][2][2] += szz * 0.5;
#endif
#endif



    }
  /* ----------------------------------------- */
  /* Assign the accumulated force from         */
  /* nonbonded interactions to the force on    */
  /* the atom.                                 */
  /* ----------------------------------------- */
    ff[k][a].x += fxa;
    ff[k][a].y += fya;
    ff[k][a].z += fza;

  }

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bonding          */
  /* contribution to the virial tensor to pvir[k].bond[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pvir[k].nbond[0] = wxx;
  pvir[k].nbond[1] = wyx;
  pvir[k].nbond[2] = wyy;
  pvir[k].nbond[3] = wzx;
  pvir[k].nbond[4] = wzy;
  pvir[k].nbond[5] = wzz;

#endif

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, assign the accumulated energies to "energy" array. */
  /*                                                                    */
  /* ================================================================== */

  energy[0] = enonbond;                 // non-shifted energy
  energy[1] = enonbond - enonbonds;     // shifted energy
  energy[2] = ecoulomb;
}

#endif
