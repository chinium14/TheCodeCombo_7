/* ======================================================================== */
/* vblock.cpp                                                               */
/*                                                                          */
/*		This subroutine uses the block accumulated values for the           */
/* properties calculated in calcvalue(), the integrates, and kinet() and    */
/* calculates the block average. It then accumulates this block average     */
/* in the resb structure to get the average of the blocks and the error.    */
/* The accumulation of the block averages is done in svalues.  svalues is   */
/* called at every iteration, while vblock is called the the frequency      */
/* sim.blockd.																*/
/*                                                                          */
/* Passed Parameters:                                                       */
/*						icyc:		The current iteration number            */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void vblock (int icyc)
{
#ifndef DOS
#ifndef MMDOS
//  double ibloc = 1.0 / sim.blockd;
  for(int k=0; k<sim.NB; k++) {

  /* ================================================================== */
  /*                                                                    */
  /* Accumulate the counter.  This value is used in outend() to get the */
  /* average of the block averages.                                     */
  /*                                                                    */
  /* ================================================================== */
	resb[k].count++;

	double icount = 1.0/(double)res[k].count;
	  
  /* ================================================================== */
  /*                                                                    */
  /* Calculate the block averages and accumulate it.                    */
  /*                                                                    */
  /* ================================================================== */
	
  /* ----------------------------------------------- */
  /* Order Parameters                                */
  /* ----------------------------------------------- */    
	ordparam[k].hela	+= ordparam[k].helb	* icount;
	ordparam[k].d_nca	+= ordparam[k].d_ncb* icount;
	ordparam[k].cona	+= ordparam[k].conb	* icount;
	ordparam[k].gyra    += ordparam[k].gyrb	* icount;
	ordparam[k].rmsda	+= ordparam[k].rmsdb* icount;
	ordparam[k].x1a		+= ordparam[k].x1b	* icount;
	ordparam[k].x2a		+= ordparam[k].x2b	* icount;

  /* ----------------------------------------------- */
  /* Temperature                                     */
  /* ----------------------------------------------- */    
	double temp  = res[k].tempb*icount;
//	double temp2 = res[k].tempc*icount;
	resb[k].tempb += temp;
//    resb[k].tempc += temp2;
        resb[k].tempc += temp*temp;


  /* ----------------------------------------------- */
  /* Pressure                                        */
  /* ----------------------------------------------- */    
#ifdef PRESSURE
    double press   = res[k].pressb*icount;
//	double press2  = res[k].pressc*icount;
	resb[k].pressb += press;
//    resb[k].pressc += press2;
    resb[k].pressc += press*press;




  /* ----------------------------------------------- */
  /* Presure tensor                                  */
  /* ----------------------------------------------- */    
    double cpress[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
//    double cpress2[6]= {0.0,0.0,0.0,0.0,0.0,0.0};
	for(int i=0; i<6; i++) {
      cpress[i] = res[k].cpressb[i]*icount;
//	  cpress2[i] = res[k].cpressc[i]*icount;
      resb[k].cpressb[i] += cpress[i];
//      resb[k].cpressc[i] += cpress2[i];
      resb[k].cpressc[i] += cpress[i]*cpress[i];
    }
#endif
	
  /* ----------------------------------------------- */
  /* Density                                         */
  /* ----------------------------------------------- */    
	double dens  = res[k].densb*icount;
//	double dens2 = res[k].densc*icount;
    resb[k].densb += dens;
//    resb[k].densc += dens2;
    resb[k].densc += dens*dens;
    
    
  /* ----------------------------------------------- */
  /* Bond Energy                                     */
  /* ----------------------------------------------- */    
    double ebond   = res[k].ebondb*icount;
    //double ebond2  = res[k].ebondc*icount;   
	resb[k].ebondb += ebond;
    //resb[k].ebondc += ebond2;
    resb[k].ebondc += ebond*ebond;
    
  /* ----------------------------------------------- */
  /* Bend Energy                                     */
  /* ----------------------------------------------- */    
    double ebend   = res[k].ebendb*icount;
    //double ebend2  = res[k].ebendc*icount;
    resb[k].ebendb += ebend;
    //resb[k].ebendc += ebend2;
    resb[k].ebendc += ebend*ebend;

  /* ----------------------------------------------- */
  /* Urey-Bradley Energy                             */
  /* ----------------------------------------------- */    
#ifndef NEUTRAL
    double eurey   = res[k].eureyb*icount;
    //double eurey2  = res[k].eureyc*icount;
    resb[k].eureyb += eurey;
    //resb[k].eureyc += eurey2;
    resb[k].eureyc += eurey*eurey;
#endif

  /* ----------------------------------------------- */
  /* Torsion Energy                                  */
  /* ----------------------------------------------- */    
    double etors   = res[k].etorsb*icount;
    //double etors2  = res[k].etorsc*icount;
    resb[k].etorsb += etors;
    //resb[k].etorsc += etors2;
    resb[k].etorsc += etors*etors;
    
  /* ----------------------------------------------- */
  /* Improper Torsion Energy                         */
  /* ----------------------------------------------- */    
    double eimpr   = res[k].eimprb*icount;
    //double eimpr2  = res[k].eimprc*icount;
    resb[k].eimprb += eimpr;
    //resb[k].eimprc += eimpr2;
    resb[k].eimprc += eimpr*eimpr;

  /* ----------------------------------------------- */
  /* Coulombic Energy                                */
  /* ----------------------------------------------- */    
	double ecoulomb  = res[k].ecoulombb*icount;
	//double ecoulomb2 = res[k].ecoulombc*icount;
	resb[k].ecoulombb += ecoulomb;
	//resb[k].ecoulombc += ecoulomb2;
	resb[k].ecoulombc += ecoulomb*ecoulomb;
	
  /* ----------------------------------------------- */
  /* Nonbonded Energy (unshifted)                    */
  /* ----------------------------------------------- */    
    double enbond   = res[k].enbondb*icount;
    //double enbond2  = res[k].enbondc*icount;
    resb[k].enbondb += enbond;
    //resb[k].enbondc += enbond2;
    resb[k].enbondc += enbond*enbond;

  /* ----------------------------------------------- */
  /* Nonbonded Energy (shifted)                      */
  /* ----------------------------------------------- */    
    double enbonds   = res[k].enbondsb*icount;
    //double enbonds2   = res[k].enbondsc*icount;
    resb[k].enbondsb += enbonds;
    //resb[k].enbondsc += enbonds2;
    resb[k].enbondsc += enbonds*enbonds;
    
  /* ----------------------------------------------- */
  /* Total Potential Energy (unshifted)              */  
  /* ----------------------------------------------- */    
    double etpot   = res[k].etpotb*icount;
    //double etpot2  = res[k].etpotc*icount;
    resb[k].etpotb += etpot;
    //resb[k].etpotc += etpot2;
    resb[k].etpotc += etpot*etpot;

  /* ----------------------------------------------- */
  /* Total Potential Energy (shifted)                */  
  /* ----------------------------------------------- */    
    double etpots   = res[k].etpotsb*icount;
    double etpots2  = res[k].etpotsc*icount; //this is the block average of PE^2
    resb[k].etpotsb += etpots;
    //resb[k].etpotsc += etpots2;
    resb[k].etpotsc += etpots*etpots;

  /* ----------------------------------------------- */
  /* Head Capacity                                   */  
  /* ----------------------------------------------- */    
    double heat_capacity = (etpots2 - etpots * etpots) / (RG*0.001) / temp / temp;
    resb[k].Cvb += heat_capacity;
    resb[k].Cvc += heat_capacity * heat_capacity;

  /* ----------------------------------------------- */
  /* Total Kinetic Energy                            */  
  /* ----------------------------------------------- */    
    double etkin   = res[k].etkinb*icount;
    //double etkin2  = res[k].etkinc*icount;
    resb[k].etkinb += etkin;
    //resb[k].etkinc += etkin2;
    resb[k].etkinc += etkin*etkin;

  /* ----------------------------------------------- */
  /* Total Energy (unshifted)                        */  
  /* ----------------------------------------------- */    
    double etotal   = res[k].etotalb*icount;
    //double etotal2  = res[k].etotalc*icount;
    resb[k].etotalb += etotal;
    //resb[k].etotalc += etotal2;
    resb[k].etotalc += etotal*etotal;

  /* ----------------------------------------------- */
  /* Total Potential Energy (shifted)                */  
  /* ----------------------------------------------- */    
    double etotals   = res[k].etotalsb*icount;
    //double etotals2  = res[k].etotalsc*icount;
    resb[k].etotalsb += etotals;
    //resb[k].etotalsc += etotals2;
    resb[k].etotalsc += etotals*etotals;

  /* ----------------------------------------------- */
  /* Extended Hamiltonian (If the simulation mode    */  
  /* is set to use a Nose Hoover thermostat or       */
  /* Parinello Rahman NPT.)                          */
  /* ----------------------------------------------- */    
	if(sim.ID == 6 || ((sim.ID2 == 4 || sim.ID2 == 5 || sim.ID2==8) && sim.ID ==1))
	{
		double H    = res[k].Hb*icount;
		double H2   = res[k].Hc*icount;
		resb[k].Hb += H;
		resb[k].Hc += H2;
	}
  }
#endif // not MMDOS
#endif //not DOS
}
