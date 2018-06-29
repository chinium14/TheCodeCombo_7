/* ======================================================================== */
/* forces.cpp                                                               */
/*                                                                          */
/*		This subroutine is the driver for the force calculations if the     */
/* multiple time step (RESPA) integration is not selected in simul.input.   */
/* It calls the subroutines that calculate the different contributions to   */
/* the forces.  It zero's out the forces on the atoms (e.g. ff[k][i].x) so  */
/* that the contributions to the different energy modes can be added in     */
/* their respective subroutines (e.g. cbond()).  It also zero's out the     */
/* virial accumulators.                                                     */
/*                                                                          */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"


  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
double cbond (int);
double cbend (int);
#ifndef NEUTRAL
double curey(int);
#endif
double ctorsion (int);
double cimproper (int);
#ifndef NLIST
void   cnonbond (int,double*);
#endif
#ifdef NLIST
void   cnonbondnl (int,double*);
#ifdef SASA								// SASA is only defined when neighbor list is defined
  double csasa(int);
#ifdef SASAREX
  double csasab(int);
#endif
#endif
#endif



#ifdef CONFIGT
double cbond2 (int);
#ifdef REST
double crestraint2(int);
#endif
double cbend2 (int);
#ifndef NEUTRAL
double curey2(int);
#endif
double ctorsion2(int);
double cimproper2(int);

#ifdef NLIST
void   cnonbondnl2 (int, double*);
#ifdef SASA
double csasa2(int);
#endif//SASA
#endif//NLIST
#endif//CONFIGT

#ifdef REST
double crestraint(int);
#endif

#ifdef SMD
double csmd(int);
#endif

#ifdef GOLIK
  #ifndef CONFIGT
    double cbb_golik(int);
    #ifdef GOBT
      double cba_golik(int);
      double cda_golik(int);
    #endif
    #ifndef DNA_GOLIK
      void cnbnd_golik(int, double*);
    #else
      void cnbnd_dnagolik(int,double*);
    #endif
  #else
    double cbb_golik2(int);
    #ifdef GOBT
      double cba_golik2(int);
      double cda_golik2(int);
    #endif
    #ifndef DNA_GOLIK
      void cnbnd_golik2(int, double*);
    #else
      void cnbnd_dnagolik2(int,double*);
    #endif
  #endif
#endif
#ifdef WALL
#ifndef CONFIGT
double cwall(int);
#else
double cwall2(int);
#endif
#endif

#ifdef CAVITY
double ccavity(int);
#endif
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

  void forces (int ibox)
{
  int k = ibox;
  /* ================================================================== */
  /*                                                                    */
  /* Zero out atomic forces arrays.                                     */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<box[k].boxns; i++) {
    ff[k][i].x = 0.0;
    ff[k][i].y = 0.0;
    ff[k][i].z = 0.0;
  }
#ifdef PR_NPT
  for(int i=0; i<9; i++) frab[k][i] = 0.0;
#endif
  /* ================================================================== */
  /*                                                                    */
  /* Zero out the virial accumulators                                   */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  for(int i=0; i<6; i++) {
    pvir[k].nbond[i] = 0.0;
    pvir[k].bond[i]  = 0.0;
    pvir[k].bend[i]  = 0.0;
#ifndef NEUTRAL
	pvir[k].urey[i]	 = 0.0;
#endif
    pvir[k].tors[i]  = 0.0;
	pvir[k].impr[i]  = 0.0;
#ifdef SASA
	pvir[k].sasa[i]  = 0.0;
#endif//SASA
#ifdef EWALD
	pvir[k].ewald_real[i] = 0.0;
	pvir[k].ewald_kspace[i] = 0.0;
	pvir[k].ewald_intra[i] = 0.0;
#endif//EWALD
  }
#endif//PRESSURE
  
  /* ================================================================== */
  /*                                                                    */
  /* Call subroutines that calculate energy and forces from the         */
  /* different contributions to the force field.                        */
  /*                                                                    */
  /* ================================================================== */

#ifndef GOLIK
#ifndef CONFIGT
/* ------------------------------------------------------------------ */
/* Regular Force Field, No CONFIGT                                    */
/* ------------------------------------------------------------------ */
  en[k].bond = cbond(k);
  #ifdef REST
  en[k].bond += crestraint(k);
  #endif

  #ifdef SMD
  en[k].bond += csmd(k);
  #endif

  en[k].bend = cbend(k);

  #ifndef NEUTRAL
  en[k].urey = curey(k);
  #endif

  en[k].tors = ctorsion(k);
  en[k].impr = cimproper(k);

  #ifdef WALL
  en[k].impr = cwall(k);
  #endif

  double interen[8];

  #ifndef NLIST
  cnonbond(k,interen);
  #endif

  #ifdef NLIST
  cnonbondnl(k,interen);
  #endif
    
  en[k].nbond  = interen[0];
  en[k].nbonds = interen[1];

  #ifdef COULOMB
  en[k].coulomb = interen[2];
  #endif

  #ifdef EWALD
  en[k].coulomb = interen[2];
  en[k].tewald  = interen[2];
  en[k].rewald  = interen[7];
  en[k].rewalds = interen[3];
  en[k].kewald  = interen[4];
  en[k].sewald  = interen[5];
  en[k].iewald  = interen[6];
  #endif

  #ifdef SASA
  en[k].esasa = csasa(k);
  #ifdef SASAREX
  en[k].esasab = csasab(k);
  #endif
  #endif
#endif // No Config T

/* ----------------------------------------------------------	*/
/* Regular Forcefield WITH CONFIGT                            */
/* ----------------------------------------------------------	*/
#ifdef CONFIGT
  config[k].hesx=0.0;
  config[k].hesy=0.0;
  config[k].hesz=0.0;
  config[k].hesr=0.0;

  en[k].bond = cbond2(k);

  #ifdef REST
  en[k].bond += crestraint2(k);
  #endif

  en[k].bend = cbend2(k);

  #ifndef NEUTRAL
  en[k].urey = curey2(k);
  #endif 

  en[k].tors = ctorsion2(k);
  en[k].impr = cimproper2(k);
  double interen[8];
  
  #ifdef NLIST
  cnonbondnl2(k,interen);
  #endif

  en[k].nbond  = interen[0];
  en[k].nbonds = interen[1];

  #ifdef COULOMB
  en[k].coulomb = interen[2];
  #endif
  #ifdef SASA
  en[k].esasa = csasa2(k);
  #endif
#endif//ConfigT

#else //GOLIK

#ifndef CONFIGT
/* ---------------------------------------------------------- */
/* GOLIK (not BGO)  Forcefield WITHOUT CONFIGT                */
/* ---------------------------------------------------------- */
  double interen[12];
  en[k].bond = cbb_golik(k);	//bond

  #ifdef GOBT
  en[k].bond += cba_golik(k);	//bend


  en[k].bond += cda_golik(k);	//tors


  #endif
  
  #ifndef DNA_GOLIK
  cnbnd_golik(k,interen);
  en[k].nbond	= interen[0];
  en[k].nbonds	= interen[1];
  en[k].bend	= interen[2];	// EPP
  en[k].urey	= interen[3];	// EPS
  en[k].impr	= interen[4];   // ESS
  #else
  cnbnd_dnagolik(k,interen);

  en[k].nbond	= interen[0];
  en[k].nbonds	= interen[1];
  en[k].bend	= interen[2];	// enative
  en[k].urey	= interen[3];	// enonnative+emm
  en[k].impr	= interen[4]; // ebp
  en[k].coulomb = interen[5];
  en[k].tors	= interen[11];
    #ifdef EWALD
    en[k].tewald  = interen[6]-interen[7]+interen[8]-interen[9]-interen[10];
    en[k].rewald  = interen[6];
    en[k].rewalds = interen[7];
    en[k].kewald  = interen[8];
    en[k].sewald  = interen[9];
    en[k].iewald  = interen[10];
    #endif
  #endif //DNA_GOLIK
  
  #ifdef WALL
  en[k].tors	= cwall(k);
  #endif

  #ifdef CAVITY
  en[k].tors	= ccavity(k);
  #endif
  #ifdef REST
  en[k].bond += crestraint(k);
  #endif

#else //CONFIGT

/* ---------------------------------------------------------- */
/* GOLIK (not BGO)  Forcefield WITH CONFIGT                   */
/* ---------------------------------------------------------- */
  config[k].hesx=0.0;
  config[k].hesy=0.0;
  config[k].hesz=0.0;
  config[k].hesr=0.0;

  double interen[6];
  en[k].bond = cbb_golik2(k);	//bond

  #ifdef GOBT
  en[k].bond += cba_golik2(k);	//bend
  en[k].bond += cda_golik2(k);	//tors
  #endif

  #ifndef DNA_GOLIK
  cnbnd_golik2(k,interen);
  en[k].nbond	= interen[0];
  en[k].nbonds	= interen[1];
  en[k].bend	= interen[2];	// EPP
  en[k].urey	= interen[3];	// EPS
  en[k].impr	= interen[4];   // ESS
  #else
  cnbnd_dnagolik2(k,interen);
  en[k].nbond	= interen[0];
  en[k].nbonds	= interen[1];
  en[k].bend	= interen[2];	// enative
  en[k].urey	= interen[3];	// enonnative+emm
  en[k].impr	= interen[4]; // ebp
  en[k].coulomb = interen[5];
  #endif

  #ifdef WALL
  en[k].tors	= cwall2(k);
  #endif

  #ifdef CAVITY
  en[k].tors	= ccavity(k);
  #endif

  #ifdef REST
  en[k].bond += crestraint2(k);
  #endif
#endif//CONFIGT
#endif //GOLIK

  /* ----------------------------------------------------------	*/
  /* If FLIM is defined, check to make sure the forces are not  */
  /* too great.  If they are, reduce them to FORCE_LIMIT.       */
  /* ----------------------------------------------------------	*/
#ifdef FLIM
  for(int i=0; i<box[k].boxns; i++){
	if(ff[k][i].x>FORCE_LIMIT) ff[k][i].x = FORCE_LIMIT;
	else if (ff[k][i].x<-FORCE_LIMIT) ff[k][i].x = -FORCE_LIMIT;
	if(ff[k][i].y>FORCE_LIMIT) ff[k][i].y = FORCE_LIMIT;
	else if (ff[k][i].y<-FORCE_LIMIT) ff[k][i].y = -FORCE_LIMIT;
	if(ff[k][i].z>FORCE_LIMIT) ff[k][i].z = FORCE_LIMIT;
	else if (ff[k][i].z<-FORCE_LIMIT) ff[k][i].z = -FORCE_LIMIT;
  }
#endif//FLIM
  
  /* ----------------------------------------------------------	*/
  /* If KONS is defined, set the forces and velocities of the  */
  /* constrained atoms to zero.                                 */
  /* ----------------------------------------------------------	*/
#ifdef KONS
  for(int i=0; i<cons[k].n; i++){
	int j = cons[k].site[i];
	ff[k][j].x = 0.0;
	ff[k][j].y = 0.0;
	ff[k][j].z = 0.0;
  }
#endif

  /* ----------------------------------------------------------	*/
  /* FORCETEST is used to compute the numerical derivatives of  */
  /* the potential energy to see if the analytical forces are   */
  /* correct.                                                   */
  /* ----------------------------------------------------------	*/
#ifdef FORCETEST
  int trial_atom = 88;
  for(int i=0; i<box[k].boxns; i++){
	  ff[k][i].x=0.0;
	  ff[k][i].y=0.0;
	  ff[k][i].z=0.0;
  }
  double energy_nb[11];
  cnbnd_dnagolik(k,energy_nb);
	//double E0 = cwall(k);;  
  
  double fxa = ff[k][trial_atom].x;
  double fya = ff[k][trial_atom].y;
  double fza = ff[k][trial_atom].z;

  double E0 = energy_nb[1];
  double E2 = energy_nb[5];
  atom[k][trial_atom].y += 0.01;

  for(int i=0; i<box[k].boxns; i++){
	  ff[k][i].x=0.0;
	  ff[k][i].y=0.0;
	  ff[k][i].z=0.0;
 }
 //double E1 = cwall(k);
  cnbnd_dnagolik(k,energy_nb);

  double fxa1 = ff[k][trial_atom].x;
  double fya1 = ff[k][trial_atom].y;
  double fza1 = ff[k][trial_atom].z;

  double fx = (fxa+fxa1)/2.0;
  double fy = (fya+fya1)/2.0;
  double fz = (fza+fza1)/2.0;
  double E1 = energy_nb[1];
  double E3 = energy_nb[5];
 
  double ftemp = -(E1-E0)/0.01;
#endif //FORCETEST
}
