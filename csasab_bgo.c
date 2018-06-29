/* ======================================================================== */
/* ======================================================================== */
/* csasab_bgo.cpp                                                            */
/*                                                                          */
/* This subroutine calculates the coarse grained (BGO) solvent energies and */
/* forces using the SASA model. It returns the potential energy due to the  */
/* solvent to be stored in en[k].esasa. It is called from either            */
/* forces() or force_short(). This subroutine is called only if SASA is     */
/* defined.  This subroutine is only called if SASA is defined.             */
/*                                                                          */
/* See: Luigi Cavallo, Jens Kleinjung and Franca Fraternali.,               */
/*      "POPS: a fast algorithm for solvent accessible			    */
/*           surface areas at atomic and residue level."                    */
/*      Nucleic Acids Research, 2003, Vol. 31, No. 13                       */
/*                                                                          */
/*      Ferrara, Philippe et. al., "Evaluation of Fast Implicit Solvent     */
/*      Model for Molecular Dynamics Simulaitons," PROTEINS: Structure      */
/*      Function, and Genetics, 46:24-33 (2002).                            */
/*                                                                          */
/* Passed Parameters:                                                       */
/*		ibox:		The box number                              */
/*                                                                          */
/* ======================================================================== */


#include "defines.h"
#ifdef SASA
#ifdef BGO
#ifdef SASAREX
#ifdef PR_NPT
void min_image_npt_full(int, double*, double*, double*);
#endif

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

double csasab (int ibox)
{
  int k	= ibox;
  int swap_diff = 0;
#ifdef SASAREX
  if (mpi.my_rank!=mpi.p-1 && mpi.flag ==1){
      swap_diff = 1;
  }
  else if (mpi.my_rank!=0 && mpi.flag==0){
      swap_diff = -1;
  }
#endif  
  
  /* ================================================================== */
  /*                                                                    */
  /* Zero out the energy accumulator                                    */
  /*                                                                    */
  /* ================================================================== */
  double esasa  = 0.0;
#ifdef WWALL
  sigW = 0.0; 
#endif
#ifndef LJ
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
  /* Get and assign information for pbc and SASA parameters.            */
  /*                                                                    */
  /* ================================================================== */
  for(int a=0; a<box[k].boxns; a++)  sasa[k][a].A =sasa[k][a].S;
  double rp   = RPROBE;
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;
 

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the energy associated with the solvent.  The first loop  */
  /* is over atoms/sites 0 to N-1 and the second loop is over its       */
  /* neighbors.                                                         */
  /*                                                                    */
  /* Note: This loop cannot be combined with the force loop because the */
  /* force needs the atomic solvation parameter which can only be       */
  /* calculated after this loop is finished.                            */
  /*                                                                    */
  /* ================================================================== */
  for(int a=0; a<box[k].boxns-1; a++) {

   
  /* ----------------------------------------------- */
  /* Assign position of particle a.                  */
  /* ----------------------------------------------- */    
    double ax = atom[k][a].x;
    double ay = atom[k][a].y;
    double az = atom[k][a].z;
   // beta[k][a] = 1.0; // initialize beta for each a. 
  /* ----------------------------------------------- */
  /* Loop around neighbors of particle a.            */
  /* ----------------------------------------------- */


    for(int bb=0; bb<nlist[k].count_sasa[a]; bb++) {
      int b = nlist[k].list_sasa[a][bb];
     
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
#ifndef PR_NPT
      if(dx >  hx) dx -= lx;
      else if(dx < -hx) dx += lx;
      if(dy >  hy) dy -= ly;
      else if(dy < -hy) dy += ly;
      if(dz >  hz) dz -= lz;
      else if(dz < -hz) dz += lz;
#endif
#ifdef PR_NPT
	min_image_npt_full(k, &dx, &dy, &dz);
#endif
  /* ----------------------------------------------- */
  /* Calculate the distance between particles a & b. */
  /* If the distance is less than the cutoff,        */
  /* accumulate the atomic solvation parameter, A.   */
  /* If the distance is greater, go to the next pair */ 
  /* in the loop.                                    */
  /* ----------------------------------------------- */
      double dr2	 = dx*dx + dy*dy + dz*dz;
	  double ra		 = sasa[k][a].r;
	  double rb		 = sasa[k][b].r;
	  double pa		 = sasa[k][a].p;
	  double pb		 = sasa[k][b].p;
	  double Sa		 = sasa[k][a].S;
	  double Sb		 = sasa[k][b].S;
	  double alpha_1 = PI*(ra+rp);			double alpha_b_1 = PI*(rb+rp);
	  double alpha_2 = (ra+rb+2.0*rp);
	  double alpha_3 = rb-ra;
	  double rc		 = alpha_2;
	  double rc2	 = rc*rc;

      if(dr2 > rc2) continue;
	  double rab = sqrt(dr2);
//        beta[k][a] += pow((rab - (ram+rbm)/(2*rp)),(1/3));
                    
#ifndef STYPE
      double pab    = ljset[k].pot[a][b].pab;
#endif

	  double bab	= alpha_1*(alpha_2- rab)*(1.0 +alpha_3/rab);
	  double bba	= alpha_b_1*(alpha_2- rab)*(1.0 -alpha_3/rab);
	  double cab	= pa*pab/Sa;
	  double cba	= pb*pab/Sb;

	  sasa[k][a].A  *= (1.0-cab*bab);
	  sasa[k][b].A  *= (1.0-cba*bba);
	}// b loop ends here

#ifdef WWALL
      double rpw = RPROBE;
      double raw = sasa[k][a].r;
      double paw = sasa[k][a].p;
      double Saw = sasa[k][a].S; 
   #ifndef STYPE
      double pabw    = 0.58720;
   #endif
      sasa[k][a].A  *= (1- paw*pabw*2*PI*(raw+rpw)*(raw-az+2*rpw)/Saw); //Spherical Crown with a water layer 
                                                    //around residue a and the on the surface

/*
      double Sbw = PI * rbw * rbw; // Flat surface area with the same radii as the residue.
      double alpha_w1 = PI*(raw+rpw);			double alpha_b_w1 = PI*(rbw+rpw);
      double alpha_w2 = (raw+rbw+2.0*rpw);
      double alpha_w3 = rbw-raw;
      double rcw		 = alpha_w2;
      double rcw2	 = rcw*rcw;
      double dxw = ax;
      double dyw = ay;
      double dzw = az + 1000;
      double drw2	 = dxw*dxw + dyw*dyw + dzw*dzw;
      if(drw2 > rcw2) continue;
      double rabw = sqrt(drw2);
   #ifndef STYPE
      double pabw    = 0.58720;
   #endif
      
      double babw	= alpha_w1*(alpha_w2- rabw)*(1.0 +alpha_w3/rabw);
      double bbaw	= alpha_b_w1*(alpha_w2- rabw)*(1.0 -alpha_w3/rabw);
      double cabw	= paw*pabw/Saw;
      double cbaw	= pbw*pabw/Sbw;

      sasa[k][a].A  *= (1.0-cabw*babw);
*/
#endif

#ifdef SPHERE
      //sasa[k][a].A *= (1.0 - )
#endif
  }// a loop ends here

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the force associated with the solvent.  The first loop   */
  /* is over atoms/sites 0 to N-1 and the second loop is over its       */
  /* neighbors.                                                         */
  /*                                                                    */
  /* ================================================================== */
  for(int a=0; a<box[k].boxns-1; a++) {
    double mBB = tfe[k][a].m_BB; 
    double mSC = tfe[k][a].m_SC; 
    double bBB = tfe[k][a].b_BB; 
    double bSC = tfe[k][a].b_SC; 
    double sBB = tfe[k][a].S_BB; 
    double sSC = tfe[k][a].S_SC; 
  // temp adjustment
  
  //  mBB = mBB * sSC / sBB;
  //  bBB = bBB * sSC / sBB;
  
  //     

  //  double alpha = sBB / sSC;
 //   double sigA = ((mBB / beta[k][a] + mSC) * concentration + (bBB / beta[k][a] + bSC)) * ((alpha + 1) * beta[k][a] / (alpha + beta[k][a])) * (sasa[k][a].A / (sBB + sSC));
#ifdef SASAREX
    double sigA = scale_factor*((mBB + mSC) * concentration[mpi.my_rank+swap_diff] + scale_factor_b*(bBB + bSC)) * (sasa[k][a].A / (sBB + sSC));
#endif

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
  /* Loop around neighbors of particle a.            */
  /* ----------------------------------------------- */
    for(int bb=0; bb<nlist[k].count_sasa[a]; bb++) {
      int b = nlist[k].list_sasa[a][bb];

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
      if(dx >  hx) dx -= lx;
      else if(dx < -hx) dx += lx;
      if(dy >  hy) dy -= ly;
      else if(dy < -hy) dy += ly;
      if(dz >  hz) dz -= lz;
      else if(dz < -hz) dz += lz;

  /* ----------------------------------------------- */
  /* Calculate the distance between particles a & b. */
  /* If the distance is less than the cutoff,        */
  /* calculate the force. If the distance is         */
  /* greater, go to the next pair in the loop.       */ 
  /* ----------------------------------------------- */
      double dr2	 = dx*dx + dy*dy + dz*dz;
	  double ra		 = sasa[k][a].r;
	  double rb		 = sasa[k][b].r;
	  double pa		 = sasa[k][a].p;
	  double pb		 = sasa[k][b].p;
	  double Sa		 = sasa[k][a].S;
	  double Sb		 = sasa[k][b].S;
	  double alpha_1 = PI*(ra+rp);			double alpha_b_1 = PI*(rb+rp);
	  double alpha_2 = (ra+rb+2.0*rp);
	  double alpha_3 = rb-ra;
	  double rc		 = alpha_2;
	  double rc2	 = rc*rc;

      if(dr2 > rc2) continue;
	  double rab = sqrt(dr2);

#ifndef STYPE
      double pab    = ljset[k].pot[a][b].pab;
#endif

	  double bab	= alpha_1	*(alpha_2- rab)*(1.0 +alpha_3/rab);
	  double bba	= alpha_b_1	*(alpha_2- rab)*(1.0 -alpha_3/rab);
	  double cab	= pa*pab/Sa;
	  double cba	= pb*pab/Sb;
	  double dbab	= -alpha_1  *(1.0 + alpha_2*alpha_3/dr2);
	  double dbba	= -alpha_b_1*(1.0 - alpha_2*alpha_3/dr2);

          double mBBb = tfe[k][b].m_BB; 
          double mSCb = tfe[k][b].m_SC; 
          double bBBb = tfe[k][b].b_BB; 
          double bSCb = tfe[k][b].b_SC; 
          double sBBb = tfe[k][b].S_BB; 
          double sSCb = tfe[k][b].S_SC; 
//          double alphab = sBBb / sSCb;
// temp ajustment           
   //       mBBb = mBBb * sSC / sBB;
   //       bBBb = bBBb * sSC / sBB;

//
          //double sigB_const = ((mBBb / beta[k][b] + mSCb) * concentration + (bBBb / beta[k][b] + bSCb)) * ((alpha + 1) * beta[k][b] / (alpha + beta[k][b])) * (1.0 / (sBBb + sSCb));
#ifdef SASAREX
          double sigB_const = scale_factor*((mBBb + mSCb) * concentration[mpi.my_rank+swap_diff] + scale_factor_b*(bBBb + bSCb)) * (1.0 / (sBBb + sSCb));
#endif
	  double force	= (sigA*(cab/(1.0-cab*bab))*dbab + 
					   sasa[k][b].A * sigB_const*(cba/(1.0-cba*bba))*dbba)/rab;

  /* ----------------------------------------- */
  /* Accumulate the forces.                    */
  /* ----------------------------------------- */
 //     fx          = force * dx;
 //     fxa        += fx;
 //    ff[k][b].x -= fx;

 //     fy          = force * dy;
 //     fya        += fy;
 //     ff[k][b].y -= fy;

 //     fz          = force * dz;
 //     fza        += fz;
 //     ff[k][b].z -= fz;

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
    }// b loop ends here

#ifdef WWALL
      double rpw = RPROBE;
      double raw = sasa[k][a].r;
      double rbw = raw;
      double paw = sasa[k][a].p;
 //     double Saw = sasa[k][a].S; 
      double dzw = az;
   #ifndef STYPE
      double pabw    = 0.58720;
   #endif
   if(dzw > (raw+2*rpw)) continue;
#ifdef SASAREX
      sigW = scale_factor*(mSCbw*concentration[mpi.my_rank+swap_diff] + scale_factor_b*bSCbw) * ((1-pbw*pabw*PI*((rbw+rpw)*(rbw+rpw)-(dzw-rpw)*(dzw-rpw))/((sBB+sSC)/4.0)));
      double forcewa = scale_factor*((mBB + mSC) * concentration[mpi.my_rank+swap_diff] + scale_factor_b*(bBB + bSC)) * 
                                    ((-(paw*pabw*2*PI*(raw+rpw))/(sBB+sSC)))+
                       scale_factor*(mSCbw * concentration[mpi.my_rank+swap_diff] + scale_factor_b*bSCbw) *
                                 (-2*pbw*pabw*PI*(abs(dzw-rpw))/((sBB+sSC)/4.0));
#endif
//      printf("a = %d, forcewa = %f,force1 = %f, \t, force2 = %f\n",a, forcewa, scale_factor*((mBB + mSC) * concentration[k] + scale_factor_b*(bBB + bSC)) *((-(paw*pabw*2*PI*(raw+rpw))/(sBB+sSC))), scale_factor*(mSCbw * concentration[k] + scale_factor_b*bSCbw)*(-2*pbw*pabw*(abs(dzw-rpw))/((sBB+sSC)/4.0)));

/*
      double pbw = 1;
      double Saw = sasa[k][a].S; 
      double Sbw = 4 * PI * rbw * rbw;
      double alpha_w1 = PI*(raw+rpw);			double alpha_b_w1 = PI*(rbw+rpw);
      double alpha_w2 = (raw+rbw+2.0*rpw);
      double alpha_w3 = rbw-raw;
      double rcw		 = alpha_w2;
      double rcw2	 = rcw*rcw;
      double dxw = ax;
      double dyw = ay;
      double dzw = az + 1000;
      double drw2	 = dxw*dxw + dyw*dyw + dzw*dzw;
      if(drw2 > rcw2) continue;
      double rabw = sqrt(drw2);
   #ifndef STYPE
      double pabw    = 0.58720;
   #endif
      
      double babw	= alpha_w1*(alpha_w2- rabw)*(1.0 +alpha_w3/rabw);
      double bbaw	= alpha_b_w1*(alpha_w2- rabw)*(1.0 -alpha_w3/rabw);
      double cabw	= paw*pabw/Saw;
      double cbaw	= pbw*pabw/Sbw;
      double dbabw	= -alpha_w1  *(1.0 + alpha_w2*alpha_w3/drw2);
      double dbbaw	= -alpha_b_w1*(1.0 - alpha_w2*alpha_w3/drw2);
      double Abw        = Sbw;
//          double alphab = sBBb / sSCb;

          //double sigB_const = ((mBBb / beta[k][b] + mSCb) * concentration + (bBBb / beta[k][b] + bSCb)) * ((alpha + 1) * beta[k][b] / (alpha + beta[k][b])) * (1.0 / (sBBb + sSCb));
          double sigB_const = scale_factor*((mBBbw + mSCbw) * concentration[k] + scale_factor_b*(bBBbw + bSCbw)) * (1.0 / (sBBbw + sSCbw));

	  double forcew	= (sigA*(cabw/(1.0-cabw*babw))*dbabw)/rabw;
// + Abw * sigB_const*(cbaw/(1.0-cbaw*bbaw))*dbbaw)/rabw;
//
*/
      double fzw          = forcewa * dzw;
      fza        += fzw;

#endif
  /* ----------------------------------------- */
  /* Assign the accumulated force from the     */
  /* implicit solvent interactions to the      */
  /* force on the atom.  Also accumulate       */
  /* the energy.                               */
  /* ----------------------------------------- */
//    ff[k][a].x += fxa;
//    ff[k][a].y += fya;
//    ff[k][a].z += fza;
    esasa      += sigA;
#ifdef WWALL
    esasa      += sigW;
#endif	


  }// aloop ends here
    double mBB = tfe[k][box[k].boxns-1].m_BB; 
    double mSC = tfe[k][box[k].boxns-1].m_SC; 
    double bBB = tfe[k][box[k].boxns-1].b_BB; 
    double bSC = tfe[k][box[k].boxns-1].b_SC; 
    double sBB = tfe[k][box[k].boxns-1].S_BB; 
    double sSC = tfe[k][box[k].boxns-1].S_SC; 
  //  double alpha = sBB / sSC;
  // temp adjustment
  
  //  mBB = mBB * sSC / sBB;
  //  bBB = bBB * sSC / sBB;
  
  //
#ifdef SASAREX
    esasa += scale_factor*((mBB + mSC) * concentration[mpi.my_rank+swap_diff] + scale_factor_b*(bBB + bSC)) * (sasa[k][box[k].boxns-1].A / (sBB + sSC));
#endif
  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bonding          */
  /* contribution to the virial tensor to pvir[k].bond[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pvir[k].sasa[0] = wxx;
  pvir[k].sasa[1] = wyx;
  pvir[k].sasa[2] = wyy;
  pvir[k].sasa[3] = wzx;
  pvir[k].sasa[4] = wzy;
  pvir[k].sasa[5] = wzz;

#endif
#endif // do all this only if LJ is not defined
return (esasa);

}
#endif
#endif
#endif
