/* ======================================================================== */
/* svalues.cpp                                                              */
/*                                                                          */
/*		This subroutine accumulates the instantaneous values for the        */
/* properties calculated in calcvalue(), the integrates, and kinet(). These */
/* accumulated values are those used to get block averages and are reset    */
/* after for each block.  The square of the value is also accumulated to    */
/* determine the error in the value.  The accumulators are found in the     */
/* res structure.                                                           */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		Simulation box number                   */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"
void end_to_end (int, int, int);
void helix(int);
void native_cont(int);
void rad_gyr(int);
#ifdef SMD
void smd_work(int);
#endif
#ifdef NSMD
void smd_work(int);
#endif
#ifdef WALL
double wall_angle(int,int);
#endif
#ifdef DNA_GOLIK
double calc_dist(int,int,int);
#endif
#ifdef TWHAM
void twham_accumulate(int);
#endif
#ifdef XWHAM
void xwham_accumulate(int);
#endif
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void svalues (int ibox)
{
#ifndef DOS
#ifndef MMDOS
#ifndef CTDOS
#ifndef XEDOS
	int k = ibox;

  /* ================================================================== */
  /*                                                                    */
  /* Accumulate the counter.  This value is used in vblock to get the   */
  /* block average of the accumulated values.                           */
  /*                                                                    */
  /* ================================================================== */
  res[k].count++;
  
  /* ================================================================== */
  /*                                                                    */
  /* Accumulate the values.                                             */
  /*                                                                    */
  /* ================================================================== */
  /* ================================================================== */
  /*                                                                    */
  /* Accumulate the order parameter of the protein molecule				*/
  /* Here average helicity and end to end distance are accumulated      */
  /*                                                                    */
  /* ================================================================== */
  end_to_end(k,SITE1,SITE2);
	helix(k);
	native_cont(k);
	rad_gyr(k);
#ifdef DNA_GOLIK
  ordparam[k].x1 = calc_dist(k,0,molec.lsite[0]-1);
#else
  ordparam[k].x1 = ordparam[k].force_1;
  ordparam[k].x2 = ordparam[k].force_2;
#endif

#ifdef WALL
	ordparam[k].x2 = wall_angle(k,wall[k].angle_site);
#endif
	ordparam[k].d_ncb += ordparam[k].d_nc;
	ordparam[k].helb  += ordparam[k].hel;
	ordparam[k].conb  += ordparam[k].con;
	ordparam[k].gyrb  += ordparam[k].gyr;
	ordparam[k].rmsdb += ordparam[k].rmsd;//Note: the average and rms of rmsd don't really mean anything here
	ordparam[k].x1b	  += ordparam[k].x1;
	ordparam[k].x2b   += ordparam[k].x2;

  /* ----------------------------------------------- */
  /* Temperature                                     */
  /* ----------------------------------------------- */    
  res[k].tempa  += box[k].temp;
  res[k].tempb  += box[k].temp; 
  res[k].tempc  += box[k].temp*box[k].temp;

  /* ----------------------------------------------- */
  /* Pressure                                        */
  /* ----------------------------------------------- */    
#ifdef PRESSURE
  res[k].pressa += box[k].press;
  res[k].pressb += box[k].press;
  res[k].pressc += box[k].press*box[k].press;

  /* ----------------------------------------------- */
  /* Presure tensor                                  */
  /* ----------------------------------------------- */    
  for(int i=0; i<6; i++) {
    res[k].cpressa[i]  += box[k].cpress[i];
    res[k].cpressb[i]  += box[k].cpress[i];
    res[k].cpressc[i]  += box[k].cpress[i]*box[k].cpress[i];
  }
#endif

  /* ----------------------------------------------- */
  /* Density                                         */
  /* ----------------------------------------------- */    
  res[k].densa		+= box[k].dens;
  res[k].densb		+= box[k].dens;
  res[k].densc		+= box[k].dens*box[k].dens;

  /* ----------------------------------------------- */
  /* Bond Energy                                     */
  /* ----------------------------------------------- */    
  res[k].ebonda		+= en[k].bond;
  res[k].ebondb		+= en[k].bond;
  res[k].ebondc		+= en[k].bond*en[k].bond;

  /* ----------------------------------------------- */
  /* Bend Energy                                     */
  /* ----------------------------------------------- */    
  res[k].ebenda		+= en[k].bend;
  res[k].ebendb		+= en[k].bend;
  res[k].ebendc		+= en[k].bend*en[k].bend;

  /* ----------------------------------------------- */
  /* Urey-Bradley Energy                             */
  /* ----------------------------------------------- */    
#ifndef NEUTRAL
  res[k].eureya		+= en[k].urey;
  res[k].eureyb		+= en[k].urey;
  res[k].eureyc		+= en[k].urey*en[k].urey;
#endif
  /* ----------------------------------------------- */
  /* Torsion Energy                                  */
  /* ----------------------------------------------- */    
  res[k].etorsa		+= en[k].tors;
  res[k].etorsb		+= en[k].tors;
  res[k].etorsc		+= en[k].tors*en[k].tors;

  /* ----------------------------------------------- */
  /* Improper torsion Energy                         */
  /* ----------------------------------------------- */    
  res[k].eimpra		+= en[k].impr;
  res[k].eimprb		+= en[k].impr;
  res[k].eimprc		+= en[k].impr*en[k].impr;

  /* ----------------------------------------------- */
  /* Coulombic Energy                                */
  /* ----------------------------------------------- */    
  res[k].ecoulomba	+= en[k].coulomb;
  res[k].ecoulombb	+= en[k].coulomb;
  res[k].ecoulombc	+= en[k].coulomb*en[k].coulomb;

  /* ----------------------------------------------- */
  /* Nonbonded Energy (unshifted)                    */
  /* ----------------------------------------------- */    
  res[k].enbonda	+= en[k].nbond;
  res[k].enbondb	+= en[k].nbond;
  res[k].enbondc	+= en[k].nbond*en[k].nbond;

  /* ----------------------------------------------- */
  /* Nonbonded Energy (shifted)                      */
  /* ----------------------------------------------- */    
  res[k].enbondsa	+= en[k].nbonds;
  res[k].enbondsb	+= en[k].nbonds;
  res[k].enbondsc	+= en[k].nbonds*en[k].nbonds;

  /* ----------------------------------------------- */
  /* Total Potential Energy (unshifted)              */  
  /* ----------------------------------------------- */    
  res[k].etpota		+= en[k].poten;
  res[k].etpotb		+= en[k].poten;
  res[k].etpotc		+= en[k].poten*en[k].poten;

  /* ----------------------------------------------- */
  /* Total Potential Energy (shifted)                */  
  /* ----------------------------------------------- */    
  res[k].etpotsa	+= en[k].potens;
  res[k].etpotsb	+= en[k].potens;
  res[k].etpotsc	+= en[k].potens*en[k].potens;

  /* ----------------------------------------------- */
  /* Total Kinetic Energy                            */  
  /* ----------------------------------------------- */    
  res[k].etkina		+= en[k].kinet;
  res[k].etkinb		+= en[k].kinet;
  res[k].etkinc		+= en[k].kinet*en[k].kinet;
 
  /* ----------------------------------------------- */
  /* Total Energy (unshifted)                        */  
  /* ----------------------------------------------- */    
  res[k].etotala	+= en[k].total;
  res[k].etotalb	+= en[k].total;
  res[k].etotalc	+= en[k].total*en[k].total;

  /* ----------------------------------------------- */
  /* Total Potential Energy (shifted)                */  
  /* ----------------------------------------------- */    
  res[k].etotalsa	+= en[k].totals;
  res[k].etotalsb	+= en[k].totals;
  res[k].etotalsc	+= en[k].totals*en[k].totals;

  /* ----------------------------------------------- */
  /* Extended Hamiltonian (If the simulation mode    */  
  /* is set to use a Nose Hoover thermostat.)        */
  /* ----------------------------------------------- */    
  if(sim.ID == 6 || ((sim.ID2 == 4 || sim.ID2 == 5 || sim.ID2==8) && sim.ID ==1))
  {
	res[k].Ha         += en[k].H;
	res[k].Hb         += en[k].H;
	res[k].Hc         += en[k].H*en[k].H;
  }
  /* The following added for sampling phi and psi	*/
  for(int i=0; i<torsN[k]; i++) {
	  tors[k][i].phia += tors[k][i].phi;
	  tors[k][i].psia += tors[k][i].psi;
	  tors[k][i].thetaa += tors[k][i].theta;
	  tors[k][i].count ++;
  }
#ifdef SMD
	smd_work(k);
#endif//SMD
#ifdef NSMD
	smd_work(k);
#endif//NSMD

  #ifdef TWHAM
  twham_accumulate(k);
  #endif
  #ifdef XWHAM
  xwham_accumulate(k);
  #endif

#endif	// NO XEDOS
#endif	// NOK CTDOS
#endif	// NOT MMDOS
#endif	// NOT WLDOS

}
