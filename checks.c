/* ======================================================================== */

/* checks.cpp                                                               */
/*                                                                          */
/*		This subroutine contains consistancy checks so that you don't       */
/* select input or define options that are not compatible.                  */
/*                                                                          */
/* ======================================================================== */
#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void checks(void)
{

  /* ================================================================== */
  /* NVE and sim.ID2                                                    */
  /* ================================================================== */
	if(sim.ID == 3){
		if(sim.ID2 > 3){
			fprintf(stdout, "If you want an NVE simulation, sim.ID2 must be 1, 2, or 3.\n");
			exit(111);
		}
	}

  /* ================================================================== */
  /* STYPE & LJ                                                         */
  /* ================================================================== */
#ifdef LJ
#ifdef STYPE
	fprintf(stdout,"Do not define STYPE if LJ is defined.\n");
	exit(111);
#endif
#endif

  /* ================================================================== */
  /* DOS & sim.ID                                                       */
  /* ================================================================== */
#ifdef DOS
	if(sim.ID!=2){
		fprintf(stdout,"You must use sim.ID==2 if DOS is defined.\n");
		exit(111);

	}
#endif

	if(sim.ID==2){
#ifndef DOS
		fprintf(stdout,"You must define DOS if sim.ID==2.\n");
		exit(111);
#endif
	}


  /* ================================================================== */
  /* MMDOS & sim.ID                                                     */
  /* ================================================================== */
#ifdef MMDOS
	if(sim.ID!=7){
		fprintf(stdout,"You must use sim.ID==7 if MMDOS is defined.\n");
		exit(111);
	}
#endif

	if(sim.ID==7){
#ifndef MMDOS
		fprintf(stdout,"You must define MMDOS if sim.ID==7.\n");
		exit(111);
#endif
	}
#ifdef MMDOS
	fprintf(stdout, "MMDOS should not be used!  It has not been updated for a long time.\n");
	exit(111);
#endif

  /* ================================================================== */
  /* CTDOS & sim.ID                                                       */
  /* ================================================================== */
#ifdef CTDOS
	if(sim.ID!=9){
		fprintf(stdout,"You must use sim.ID==9 if CTDOS is defined.\n");
		exit(111);
	}
#endif

	if(sim.ID==9){
#ifndef CTDOS
		fprintf(stdout,"You must define CTDOS if sim.ID==9.\n");
		exit(111);
#endif
	}


  /* ================================================================== */
  /* XEDOS & sim.ID                                                       */
  /* ================================================================== */
#ifdef XEDOS
	if(sim.ID!=10){
		fprintf(stdout,"You must use sim.ID==10 if XEDOS is defined.\n");
		exit(111);
	}
#endif

	if(sim.ID==10){
#ifndef XEDOS
		fprintf(stdout,"You must define XEDOS if sim.ID==10.\n");
		exit(111);
#endif
	}
  /* ================================================================== */
  /* FX_EDOS & sim.ID                                                       */
  /* ================================================================== */
#ifdef FX_EDOS
	if(sim.ID!=11){
		fprintf(stdout,"You must use sim.ID==11 if FX_EDOS is defined.\n");
		exit(111);

	}
#endif

	if(sim.ID==11){
#ifndef FX_EDOS
		fprintf(stdout,"You must define FX_EDOS if sim.ID==11.\n");
		exit(111);
#endif
	}
  /* ================================================================== */
  /* TDXEDOS & sim.ID                                                   */
  /* ================================================================== */
#ifdef TDXEDOS
	if(sim.ID!=12){
		fprintf(stdout,"You must use sim.ID==12 if TDXEDOS is defined.\n");
		exit(111);
	}
#endif

	if(sim.ID==12){
#ifndef TDXEDOS
		fprintf(stdout,"You must define TDXEDOS if sim.ID==12.\n");
		exit(111);
#endif
	}
  /* ================================================================== */
  /* NMA & sim.ID                                                       */
  /* ================================================================== */
#ifdef NMA
	if(sim.ID!=8){
		fprintf(stdout,"You must use sim.ID==8 if NMA is defined.\n");
		exit(111);
	}
#endif 

	if(sim.ID==8){
#ifndef NMA
		fprintf(stdout,"You must define NMA if sim.ID==8.\n");
		exit(111);
#endif
	}


  /* ================================================================== */
  /* PR_NPT & sim.ID                                                    */
  /* ================================================================== */
#ifdef PR_NPT
//	if(sim.ID!=6 || ){
//		fprintf(stdout,"You must use sim.ID==6 if PR_NPT is defined.\n");
//		exit(111);
//	}
#endif 

	if(sim.ID==6 && (sim.ID2 !=9 && sim.ID2 !=10 && sim.ID2 !=1 && sim.ID2 !=2 && sim.ID2 !=3)){
#ifndef PR_NPT
		fprintf(stdout,"You must define PR_NPT if sim.ID==6.");
		exit(111);
#endif
		if(sim.ID2!=6 && sim.ID2!=7){
			fprintf(stdout,"sim.ID2 must = 6 or 7 if sim.ID==6.\n");
			exit(111);
		}
	}

  /* ================================================================== */
  /* PR_NPT & PRESSURE                                                  */
  /* ================================================================== */
#ifdef PR_NPT
#ifndef PRESSURE
		fprintf(stdout,"You must define PRESSURE if PR_NPT if defined.\n");
		exit(111);
#endif
#endif

  /* ================================================================== */
  /* SASA & NEUTRAL, RDIE, AND COULOMB                                  */
  /* ================================================================== */
#ifdef SASA
int flag = 0;
#if !defined(NEUTRAL) && !defined(BGO)
	fprintf(stdout,"You must define NEUTRAL if SASA is defined.\n");
	flag = 1;
#endif
#if !defined(RDIE) && !defined(BGO)
	fprintf(stdout,"You must define RDIE if SASA is defined.\n");
	flag = 1;
#endif
#if !defined(COULOMB) && !defined(BGO)
	fprintf(stdout,"You must define COULOMB if SASA is defined.\n");
	flag = 1;
#endif
	if(flag == 1) exit(111);
#endif
	
  /* ================================================================== */
  /* COULOMB & EWALD                                                    */
  /* ================================================================== */
#ifdef EWALD
#ifdef COULOMB
	fprintf(stdout,"You cannot define COULOMB if EWALD is defined.\n");
	exit(111);
#endif
#endif

#ifdef COULOMB
#ifdef EWALD
	fprintf(stdout,"You cannot define EWALD if COULOMB is defined.\n");
	exit(111);
#endif
#endif

  /* ================================================================== */
  /* NLIST                                                             */
  /* ================================================================== */
#ifndef NLIST
	fprintf(stdout,"You must define NLIST or non-bonded interactions will not be calculated properly.\n");
	exit(111);
#endif
	
  /* ================================================================== */
  /* ELASTIC                                                            */
  /* ================================================================== */
#ifdef ELASTIC
#ifndef PRESSURE
	fprintf(stdout,"You must define PRESSURE if ELASTIC is defined.\n");
	exit(111);
#endif
#endif

  /* ================================================================== */
  /* PR_NPT and CLIST                                                   */
  /* ================================================================== */
#ifdef PR_NPT
#ifdef CLIST
	fprintf(stdout,"You cannot use CLIST if PR_NPT is defined.\n");
	exit(111);
#endif
#endif
#ifdef CLIST
#ifdef PR_NPT
	fprintf(stdout,"You cannot use PR_NPT if CLIST is defined.\n");
	exit(111);
#endif
#endif

  /* ================================================================== */
  /* PR_NPT and CLIST                                                   */
  /* ================================================================== */
#ifdef PIVOT
#ifdef CLIST
	fprintf(stdout,"You cannot use CLIST if PIVOT is defined.\n");
	exit(111);
#endif
#endif
#ifdef CLIST
#ifdef PIVOT
	fprintf(stdout,"You cannot use PIVOT if CLIST is defined.\n");
	exit(111);
#endif
#endif



  /* ================================================================== */
  /* CLIST and NPT                                                      */
  /* ================================================================== */
if(sim.ID==6){
#ifdef CLIST
        fprintf(stdout,"NPT is not coded to work with CLIST.\n");
        exit(111);

#endif
}
#ifdef CLIST
if(sim.ID==6){
        fprintf(stdout,"NPT is not coded to work with CLIST.\n");
        exit(111);
}
#endif

  /* ================================================================== */
  /* NPT and integration methods.                                       */
  /* ================================================================== */
if(sim.ID==6){
  if(sim.ID2 == 4 || sim.ID2 == 5 || sim.ID2 == 8){
    fprintf(stdout,"You must use either 1, 2, 3, 6, 7, 9, or 10 for sim.ID2 if sim.ID = 6\n");
    exit(111);
  }
}
  /* ================================================================== */
  /* KONS and sim.ID                                                    */
  /* ================================================================== */
#ifdef KONS
	if((sim.ID != 1) && (sim.ID != 3)){
		fprintf(stdout, "KONS has only been tested for NVE and NVT simulations. \n");
		exit(111);
	}
#endif

  /* ================================================================== */
  /* CONFIGT and REST                                                   */
  /* ================================================================== */
/*#ifdef REST
#ifdef CONFIGT
	fprintf(stdout,"You cannot use CONFIGT if REST is defined(there is no crestraint2.cpp).\n");
	exit(111);
#endif
#endif
#ifdef CONFIGT
#ifdef REST
	fprintf(stdout,"You cannot use REST if CONFIGT is defined(there is no crestraint2.cpp).\n");
	exit(111);
#endif
#endif*/

  /* ================================================================== */
  /* CONFIGT and Multiple Times Step                                    */
  /* ================================================================== */

#ifdef CONFIGT
  if(sim.ID2 == 3 || sim.ID2 == 5){
	  fprintf(stdout,"You cannot use CONFIGT with multiple time step.\n");
	  exit(111);
  }
#endif


  /* ================================================================== */
  /* DOS probability of moves check                                     */
  /* ================================================================== */
#ifdef DOS
	for(int k=0; k<sim.NB; k++){
		if((mc_pivot.PMS[k] + mc_trans.PMS[k] + mc_axial.PMS[k] + mc_rand.PMS[k]) > 1.0){
			fprintf(stdout,"The sum of the DOS probabilites (%f + %f + %f + %f = %f) in box %d is greater than 1.0.\n", mc_pivot.PMS[k], mc_trans.PMS[k], mc_axial.PMS[k], mc_rand.PMS[k], mc_pivot.PMS[k] + mc_trans.PMS[k] + mc_axial.PMS[k] + mc_rand.PMS[k], k);
			exit(111);
		}
	}
#endif
#ifdef CTDOS
	for(int k=0; k<sim.NB; k++){
		if((mc_pivot.PMS[k] + mc_trans.PMS[k] + mc_axial.PMS[k] + mc_rand.PMS[k]) > 1.0){
			fprintf(stdout,"The sum of the CTDOS probabilites (%f + %f + %f + %f = %f) in box %d is greater than 1.0.\n", mc_pivot.PMS[k], mc_trans.PMS[k], mc_axial.PMS[k], mc_rand.PMS[k], mc_pivot.PMS[k] + mc_trans.PMS[k] + mc_axial.PMS[k] + mc_rand.PMS[k], k);
			exit(111);
		}
	}
#endif
#ifdef XEDOS
	for(int k=0; k<sim.NB; k++){
        float P = mc_pivot.PMS[k] + mc_trans.PMS[k] + mc_axial.PMS[k] + mc_rand.PMS[k];
		if( P > 1.0){
			fprintf(stdout,"The sum of the XEDOS probabilites (%f + %f + %f + %f = %f) in box %d is greater than 1.0.\n", mc_pivot.PMS[k], mc_trans.PMS[k], mc_axial.PMS[k], mc_rand.PMS[k], P,k);
			exit(111);
		}
	}
#endif
  /* ================================================================== */
  /* ROTN and PR_NPT                                                    */
  /* ================================================================== */
#ifdef ROTN
#ifdef PR_NPT
	fprintf(stdout,"You cannot use PR_NPT if ROTN is defined.\n");
	exit(111);
#endif
#endif
#ifdef PR_NPT
#ifdef ROTN
	fprintf(stdout,"You cannot use ROTN if PR_NPT is defined.\n");
	exit(111);
#endif
#endif


  /* ================================================================== */
  /* MC, PRESSURE, SMD, EWALD, CTDOS                                    */
  /* ================================================================== */
#ifdef MC
//#ifdef CTDOS
//    fprintf(stdout,"You cannot use MC if CTDOS is defined.\n");
//   exit(111);
//#endif
//#ifdef PRESSURE
//	fprintf(stdout,"You cannot use PRESSURE if MC is defined.\n");
//	exit(111);
//#endif

//#ifdef EWALD
//	fprintf(stdout,"You cannot use EWALD if MC is defined.\n");
//	exit(111);
//#endif

#ifdef SMD
	fprintf(stdout,"You cannot use SMD if MC is defined.\n");
	exit(111);
#endif
#ifdef SLIST
	for(int k=0; k<sim.NB; k++){
		if (mc_solv.delta[k] > (sim.rclist - sim.rc)*0.5){
			fprintf(stdout,"You cannot use a delta (%f) greater than half the shell radius (%f).\n", mc_solv.delta[k], (sim.rclist - sim.rc)*0.5);
			exit(111);
		}
	}
#endif

#endif// MC

  /* ================================================================== */
  /* NLIST                                                              */
  /* ================================================================== */
#ifndef NLIST
	fprintf(stdout, "You must define NLIST.\n");
	exit(111);
#endif

  /* ================================================================== */
  /* BEAD and S14                                                       */
  /* ================================================================== */
#ifdef BEAD
#ifdef S14
        fprintf(stdout,"You cannot use S14 if BEAD is defined.\n");
        exit(111);
#endif
#endif
#ifdef S14
#ifdef BEAD
        fprintf(stdout,"You cannot use BEAD if S14 is defined.\n");
        exit(111);
#endif
#endif
  /* ================================================================== */
  /* GOLIK & MTS                                                        */
  /* ================================================================== */
#ifdef GOLIK
	if(sim.ID2 == 3 || sim.ID2 == 5){
        fprintf(stdout,"You cannot use GOLIK with multiple time step.\n  Change force_short and force_long if you want to do so.\n");
        exit(111);
	}
#endif

  /* ================================================================== */
  /* GOLIK and PR_NPT                                                   */
  /* ================================================================== */
#ifdef GOLIK
#ifdef PR_NPT
        fprintf(stdout,"You cannot use PR_NPT if GOLIK is defined.\n");
        exit(111);
#endif
#endif
#ifdef PR_NPT
#ifdef GOLIK
        fprintf(stdout,"You cannot use GOLIK if PR_NPT is defined.\n");
        exit(111);
#endif
#endif

  /* ================================================================== */
  /* DNA_GOLIK and STYPE                                                */
  /* ================================================================== */
#ifdef DNA_GOLIK
#ifndef STYPE
        fprintf(stdout,"You must use STYPE if DNA_GOLIK is defined.\n");
        exit(111);
#endif
#endif
  /* ================================================================== */
  /* GOLIK and STYPE                                                    */
  /* ================================================================== */
#ifdef GOLIK
#ifndef DNA_GOLIK
#ifdef STYPE
        fprintf(stdout,"You cannot use STYPE if GOLIK is defined.\n");
        exit(111);
#endif
#endif
#endif

#ifdef STYPE
#ifndef DNA_GOLIK
#ifdef GOLIK
        fprintf(stdout,"You cannot use GOLIK if STYPE is defined.\n");
        exit(111);
#endif
#endif
#endif

  /* ================================================================== */
  /* GOLIK & MTS													    */
  /* ================================================================== */
#ifdef GOLIK

		if(sim.nsteps != 1){
			fprintf(stdout,"nsteps is %i. You cannot use nsteps > 1 if GOLIK is defined.\n", sim.nsteps);
			exit(111);
		}

#endif


  /* ================================================================== */
  /* GOLIK and DNA_GOLIK                                                */
  /* ================================================================== */
#ifdef DNA_GOLIK
#ifndef GOLIK
        fprintf(stdout,"You must use GOLIK if DNA_GOLIK is defined.\n");
        exit(111);
#endif
#endif

  /* ================================================================== */
  /* GOLIK and GOBT                                                     */
  /* ================================================================== */
#ifdef GOBT
#ifndef GOLIK
        fprintf(stdout,"You must use GOLIK if GOBT is defined.\n");
        exit(111);
#endif
#endif

  /* ================================================================== */
  /* WALL and TLATE                                                     */
  /* ================================================================== */
#ifdef WALL
#ifdef TLATE
        fprintf(stdout,"You cannot use TLATE if WALL is defined.\n");
        exit(111);
#endif
#endif
#ifdef TLATE
#ifdef WALL
        fprintf(stdout,"You cannot use WALL if TLATE is defined.\n");
        exit(111);
#endif
#endif

  /* ================================================================== */
  /* ACAVITY AND BCAVITY and TLATE                                      */
  /* ================================================================== */
#if defined(ACAVITY) || defined(BCAVITY)
#ifdef TLATE
        fprintf(stdout,"You cannot use TLATE if ACAVITY or BCAVITY are defined.\n");
        exit(111);
#endif
#endif

  /* ================================================================== */
  /* WALL and MTS                                                       */
  /* ================================================================== */
#ifdef WALL
	if(sim.ID2 == 3 || sim.ID2 == 5){
		fprintf(stdout,"You cannot use multiple time step if WALL is defined.\n");
        exit(111);
	}
#endif

  /* ================================================================== */
  /* WALL and GOLIK                                                     */
  /* ================================================================== */
//#ifdef WALL
//#ifndef GOLIK
//        fprintf(stdout,"You must use GOLIK if WALL is defined.\n");
//        exit(111);
//#endif
//#endif

  /* ================================================================== */
  /* WALL and AWALL, RWALL, and NWALL                                   */
  /* ================================================================== */
#ifdef WALL
int aflag = 0;
int rflag = 0;
int nflag = 0;
int hflag = 0;
int cflag = 0;
int wflag = 0;
int sflag = 0;
#ifdef AWALL
aflag = 1;
#endif
#ifdef RWALL
rflag = 1;
#endif
#ifdef NWALL
nflag = 1;
#endif
#ifdef HWALL
hflag =1;
#endif
#ifdef CWALL
cflag =1;
#endif
#ifdef WWALL
wflag=1;
#endif
#ifdef SPHERE
sflag=1;
#endif
int sum_flag = aflag + rflag + nflag + hflag + cflag + wflag + sflag;
if(sum_flag != 1){
        fprintf(stdout,"If WALL is defined, you must also define one and only one of the following: AWALL, RWALL, NWALL, HWALL, CWALL, WWALL, SPHERE.\n");
        exit(111);

}
#elif defined(AWALL) || defined(RWALL) || defined(NWALL) || defined(HWALL) || defined(CWALL) || defined(WWALL) || defined(SPHERE)
  // If WALL was not defined, define it here
  
//  #define WALL   // <-- This would be nice, but we already read in all the input stuff

	fprintf(stdout,"You must define WALL if one of [ARNHC]WALL is defined!\n");
	exit(111);

#endif

  /* ================================================================== */
  /* WALL and AWALL, RWALL, and NWALL                                   */
  /* ================================================================== */
#ifdef WALL
#ifdef CONFIGT
  fprintf(stdout,"You can't use CONFIGT with WALL.  Need to code change cwall2 to use regular force field\nSee cwall for example.\n");
  exit(111);
#endif

#endif

  /* ================================================================== */
  /* PREEQUIL and modes                                                 */
  /* ================================================================== */
#ifdef PREEQUIL
  if(sim.ID != 1 && sim.ID != 3 && sim.ID != 5 && sim.ID != 6){
    fprintf(stdout,"You can only use sim.ID's of 1, 3, 5, or 6 if PREEQUIL is defined.\n");
    exit(111);
  }
#endif


  /* ================================================================== */
  /* RFC and CONFIGT                                                    */
  /* ================================================================== */
#ifndef DNA_GOLIK
#ifdef RFC
#ifdef CONFIGT
        fprintf(stdout,"You cannot use CONFIGT if RFC is defined.\n");
        exit(111);
#endif
#endif
#ifdef CONFIGT
#ifdef RFC
        fprintf(stdout,"You cannot use RFC if CONFIGT is defined.\n");
        exit(111);
#endif
#endif
#endif

  /* ================================================================== */
  /* SHIFTV and CONFIGT                                                 */
  /* ================================================================== */
#ifdef SHIFTV
#ifdef CONFIGT
        fprintf(stdout,"You cannot use CONFIGT if SHIFTV is defined.\n");
        exit(111);
#endif
#endif
#ifdef CONFIGT
#ifdef SHIFTV
        fprintf(stdout,"You cannot use SHIFTV if CONFIGT is defined.\n");
        exit(111);
#endif
#endif

  /* ================================================================== */
  /* DNA_GOLIK and CONFIGT                                              */
  /* ================================================================== */
#ifdef DNA_GOLIK
#ifdef CONFIGT
        fprintf(stdout,"You cannot use CONFIGT if DNA_GOLIK is defined.\n");
        fprintf(stdout,"The code for cgolik2 must be updated to include the counterions.\n");
        exit(111);
#endif
#endif

  /* ================================================================== */
  /* DNA_GOLIK and nvt_mc                                               */
  /* ================================================================== */
#ifdef DNA_GOLIK
  if(sim.ID == 4){
        fprintf(stdout,"You cannot use NVT MC (sim.ID=4) if DNA_GOLIK is defined.\n");
        fprintf(stdout,"The code for cgolik_mc, cgolik_mc2 must be updated to include the counterions.\n");
        exit(111);
  }
#endif

  /* ================================================================== */
  /* BGO and other define options.                                   */
  /* ================================================================== */
/*#ifdef BGO
#if defined(COULOMB) || defined(GOLIK) || defined(EWALD)
   fprintf(stdout,"You cannot use BGO with COULOMB, EWALD, or GOLIK.\n");
   exit(111);
#endif
#endif
*/
  /* ============================================================ */
  /* ZHOU and PRESSURE                                            */
  /* ============================================================ */
#ifdef ZHOU
#ifndef PRESSURE
fprintf(stdout, "Warning: ZHOU has no effect if PRESSURE is not defined!!!\n");
#endif
#endif


  /* ============================================================ */
  /* KONS and ZEROAM                                              */
  /* ============================================================ */
#ifdef KONS
#ifdef ZEROAM
fprintf(stdout, "Warning: KONS cannot be used if ZEROAM is defined!!!\n");
#endif
#endif

  /* ============================================================ */
  /* KONS and ZEROAM                                              */
  /* ============================================================ */
#ifdef TWHAM
if( !(sim.ID==1 || sim.ID==5) )
{
	fprintf(stdout, "ERROR: TWHAM only works with sim.ID = 1 (NVT) or 5 (replica exchange) !!!\n");
	exit(111);
}
#endif


}
