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
void dos_vblock (unsigned long icyc)
{
	/*------------------------------------------------- */
	/* If DOS is defined just collect the thermodynamic */
	/* values for each energy bin and do not accumulate */
	/* separately as above.								*/
	/*--------------------------------------------------*/
#ifdef DOS
  for(int k=0; k<sim.NB; k++) {
//	  double helsum = 0.0;
//	  double dncsum = 0.0;
//	  double consum = 0.0;
	  for(int i=0; i<sim_dos[k].e_bins; i++){
		  average[k][i].ebonda		= average[k][i].ebond	/ average[k][i].count;
		  average[k][i].ebenda		= average[k][i].ebend	/ average[k][i].count;
#ifndef NEUTRAL
		  average[k][i].eureya      = average[k][i].eurey   / average[k][i].count;
#endif
		  average[k][i].etorsa		= average[k][i].etors	/ average[k][i].count;
		  average[k][i].eimpra		= average[k][i].eimpr	/ average[k][i].count;
		  average[k][i].elja		= average[k][i].elj		/ average[k][i].count;
		  average[k][i].ecoula		= average[k][i].ecoul	/ average[k][i].count;
		  average[k][i].epotensa	= average[k][i].epotens	/ average[k][i].count;
		  average[k][i].d_nca		= average[k][i].d_nc	/ average[k][i].count;
		  average[k][i].hela		= average[k][i].hel		/ average[k][i].count;
		  average[k][i].cona		= average[k][i].con		/ average[k][i].count;
		  average[k][i].con_2a		= average[k][i].con_2	/ average[k][i].count;
		  average[k][i].gyra		= average[k][i].gyr		/ average[k][i].count;
		  average[k][i].rmsda		= average[k][i].rmsd    / average[k][i].count;
#ifdef SASA
		  average[k][i].esasaa		= average[k][i].esasa	/ average[k][i].count;
#endif

	  }

  }
#endif //DOS

  	/*--------------------------------------------------- */
	/* If MMDOS is defined just collect the thermodynamic */
	/* values for each energy bin and do not accumulate   */
	/* separately as above.								  */
	/*----------------------------------------------------*/
#ifdef MMDOS
  for(int k=0; k<sim.NB; k++) {
//	  double helsum = 0.0;
//	  double dncsum = 0.0;
//	  double consum = 0.0;
	  for(int i=0; i<sim_dos[k].e_bins; i++){
		  average[k][i].ebonda		= average[k][i].ebond	/ average[k][i].count;
		  average[k][i].ebenda		= average[k][i].ebend	/ average[k][i].count;
#ifndef NEUTRAL
		  average[k][i].eureya      = average[k][i].eurey   / average[k][i].count;
#endif
		  average[k][i].etorsa		= average[k][i].etors	/ average[k][i].count;
		  average[k][i].eimpra		= average[k][i].eimpr	/ average[k][i].count;
		  average[k][i].elja		= average[k][i].elj		/ average[k][i].count;
		  average[k][i].ecoula		= average[k][i].ecoul	/ average[k][i].count;
		  average[k][i].epotensa	= average[k][i].epotens	/ average[k][i].count;
		  average[k][i].d_nca		= average[k][i].d_nc	/ average[k][i].count;
		  average[k][i].hela		= average[k][i].hel		/ average[k][i].count;
		  average[k][i].cona		= average[k][i].con		/ average[k][i].count;
		  average[k][i].con_2a		= average[k][i].con_2	/ average[k][i].count;
		  average[k][i].rmsda		= average[k][i].rmsd    / average[k][i].count;

#ifdef SASA
		  average[k][i].esasaa		= average[k][i].esasa	/ average[k][i].count;
#endif
	  }

 }
#endif //MMDOS
#ifdef CTDOS
  for(int k=0; k<sim.NB; k++) {
//	  double helsum = 0.0;
//	  double dncsum = 0.0;
//	  double consum = 0.0;
	  for(int i=0; i<sim_dos[k].e_bins; i++){
		  average[k][i].ebonda		= average[k][i].ebond	/ average[k][i].count;
		  average[k][i].ebenda		= average[k][i].ebend	/ average[k][i].count;
#ifndef NEUTRAL
		  average[k][i].eureya      = average[k][i].eurey   / average[k][i].count;
#endif
		  average[k][i].etorsa		= average[k][i].etors	/ average[k][i].count;
		  average[k][i].eimpra		= average[k][i].eimpr	/ average[k][i].count;
		  average[k][i].elja		= average[k][i].elj		/ average[k][i].count;
		  average[k][i].ecoula		= average[k][i].ecoul	/ average[k][i].count;
		  average[k][i].epotensa	= average[k][i].epotens	/ average[k][i].count;
		  average[k][i].d_nca		= average[k][i].d_nc	/ average[k][i].count;
		  average[k][i].hela		= average[k][i].hel		/ average[k][i].count;
		  average[k][i].cona		= average[k][i].con		/ average[k][i].count;
		  average[k][i].con_2a		= average[k][i].con_2	/ average[k][i].count;
		  average[k][i].gyra		= average[k][i].gyr		/ average[k][i].count;
		  average[k][i].rmsda		= average[k][i].rmsd    / average[k][i].count;
#ifdef SASA
		  average[k][i].esasaa		= average[k][i].esasa	/ average[k][i].count;
#endif
	  }

  }
#endif //CTDOS

#ifdef XEDOS
  for(int k=0; k<sim.NB; k++) {
//	  double helsum = 0.0;
//	  double dncsum = 0.0;
//	  double consum = 0.0;
	  for(int i=0; i<sim_dos[k].l_bins; i++){
		  average[k][i].ebonda		= average[k][i].ebond	/ average[k][i].count;
		  average[k][i].ebenda		= average[k][i].ebend	/ average[k][i].count;
#ifndef NEUTRAL
		  average[k][i].eureya      = average[k][i].eurey   / average[k][i].count;
#endif
		  average[k][i].etorsa		= average[k][i].etors	/ average[k][i].count;
		  average[k][i].eimpra		= average[k][i].eimpr	/ average[k][i].count;
		  average[k][i].elja		= average[k][i].elj		/ average[k][i].count;
		  average[k][i].ecoula		= average[k][i].ecoul	/ average[k][i].count;
		  average[k][i].epotensa	= average[k][i].epotens	/ average[k][i].count;
		  average[k][i].d_nca		= average[k][i].d_nc	/ average[k][i].count;
		  average[k][i].hela		= average[k][i].hel		/ average[k][i].count;
		  average[k][i].cona		= average[k][i].con		/ average[k][i].count;
		  average[k][i].con_2a		= average[k][i].con_2	/ average[k][i].count;
		  average[k][i].force_1a	= average[k][i].force_1	/ average[k][i].count;
		  average[k][i].force_2a	= average[k][i].force_2	/ average[k][i].count;
		  average[k][i].gyra		= average[k][i].gyr		/ average[k][i].count;
		  average[k][i].rmsda		= average[k][i].rmsd    / average[k][i].count;
		  
#ifdef SASA
		  average[k][i].esasaa		= average[k][i].esasa	/ average[k][i].count;
#endif
//		  helsum	+=	average[k][i].hela;
//		  dncsum	+=	average[k][i].d_nca;
//		  consum    +=  average[k][i].cona;
	  }
//	  ordparam[k].hela  = helsum/sim_dos[k].l_bins;
//	  ordparam[k].d_nca = dncsum/sim_dos[k].l_bins;
//	  ordparam[k].cona  = consum/sim_dos[k].l_bins;
  }
#endif //XEDOS
 #ifdef FX_EDOS
  for(int k=0; k<sim.NB; k++) {
//	  double helsum = 0.0;
//	  double dncsum = 0.0;
//	  double consum = 0.0;
	  for(int i=0; i<sim_dos[k].e_bins; i++){
		  average[k][i].ebonda		= average[k][i].ebond	/ average[k][i].count;
		  average[k][i].ebenda		= average[k][i].ebend	/ average[k][i].count;
#ifndef NEUTRAL
		  average[k][i].eureya      = average[k][i].eurey   / average[k][i].count;
#endif
		  average[k][i].etorsa		= average[k][i].etors	/ average[k][i].count;
		  average[k][i].eimpra		= average[k][i].eimpr	/ average[k][i].count;
		  average[k][i].elja		= average[k][i].elj		/ average[k][i].count;
		  average[k][i].ecoula		= average[k][i].ecoul	/ average[k][i].count;
		  average[k][i].epotensa	= average[k][i].epotens	/ average[k][i].count;
		  average[k][i].d_nca		= average[k][i].d_nc	/ average[k][i].count;
		  average[k][i].hela		= average[k][i].hel		/ average[k][i].count;
		  average[k][i].cona		= average[k][i].con		/ average[k][i].count;
		  average[k][i].con_2a		= average[k][i].con_2	/ average[k][i].count;
		  average[k][i].gyra		= average[k][i].gyr		/ average[k][i].count;
		  average[k][i].rmsda		= average[k][i].rmsd    / average[k][i].count;
#ifdef SASA
		  average[k][i].esasaa		= average[k][i].esasa	/ average[k][i].count;
#endif

	  }

  }
#endif //FX_EDOS

#ifdef TDXEDOS
  for(int k=0; k<sim.NB; k++) {
	  for(int i=0; i<sim_dos[k].x1_bins; i++){
		  for(int j=0; j<sim_dos[k].x2_bins; j++){	
			average[k][i][j].ebonda		= average[k][i][j].ebond	/ average[k][i][j].count;
			average[k][i][j].ebenda		= average[k][i][j].ebend	/ average[k][i][j].count;
			#ifndef NEUTRAL
			average[k][i][j].eureya      = average[k][i][j].eurey   / average[k][i][j].count;
			#endif
			average[k][i][j].etorsa		= average[k][i][j].etors	/ average[k][i][j].count;
			average[k][i][j].eimpra		= average[k][i][j].eimpr	/ average[k][i][j].count;
			average[k][i][j].elja		= average[k][i][j].elj		/ average[k][i][j].count;
			average[k][i][j].ecoula		= average[k][i][j].ecoul	/ average[k][i][j].count;
			average[k][i][j].epotensa	= average[k][i][j].epotens	/ average[k][i][j].count;
			average[k][i][j].d_nca		= average[k][i][j].d_nc		/ average[k][i][j].count;
			average[k][i][j].hela		= average[k][i][j].hel		/ average[k][i][j].count;
			average[k][i][j].cona		= average[k][i][j].con		/ average[k][i][j].count;
			average[k][i][j].force_1a	= average[k][i][j].force_1	/ average[k][i][j].count;
			average[k][i][j].force_2a	= average[k][i][j].force_2	/ average[k][i][j].count;
			average[k][i][j].gyra		= average[k][i][j].gyr		/ average[k][i][j].count;
			average[k][i][j].rmsda		= average[k][i][j].rmsd		/ average[k][i][j].count;
			#ifdef SASA
			average[k][i][j].esasaa		= average[k][i][j].esasa	/ average[k][i][j].count;
			#endif
			average[k][i][j].x1a    = average[k][i][j].x1   / average[k][i][j].count;
			average[k][i][j].x2a    = average[k][i][j].x2   / average[k][i][j].count;
		  }
	  }
  }
#endif //TDXEDOS
}
