/* ======================================================================== */
/* dos_svalues.cpp                                                          */
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
void configtemp(int);
void native_cont(int);
void rad_gyr(int);
#ifdef WALL
double wall_angle(int, int);
#endif
#ifdef RMSD
int  quatfit			(int, int, int);
#endif
#ifdef DNA_GOLIK
double calc_dist(int,int,int);
#endif
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void dos_svalues (int ibox)
{
  int k = ibox;

		/* ------------------------------------------------------------------------------------ */
		/* Call order parameter subroutines.                                                    */
        /* ------------------------------------------------------------------------------------ */

#ifdef RMSD
  quatfit(k,rmsd_modes.file,rmsd_modes.atoms);
#endif
  helix(k);
  native_cont(k);
#ifdef TDXEDOS
  ordparam[k].x1 = ordparam[k].d_nc;
  ordparam[k].x2 = wall_angle(k,wall[k].angle_site);
#else
#ifdef WALL
  ordparam[k].x2 = wall_angle(k,wall[k].angle_site);
#endif
#ifdef DNA_GOLIK
  ordparam[k].x1 = calc_dist(k,0,molec.lsite[0]-1);
#endif
#endif  
  /* ----------	*/
  /*    DOS		*/
  /* ----------	*/
#ifdef DOS
  int bin =  (int) ((en[k].potens - sim_dos[k].e_begin)/sim_dos[k].e_width);	
  average[k][bin].count ++;
  average[k][bin].ebond		+=	en[k].bond;
  average[k][bin].ebend		+=	en[k].bend;
#ifndef NEUTRAL
  average[k][bin].eurey		+=	en[k].urey;
#endif
  average[k][bin].etors		+=	en[k].tors;
  average[k][bin].eimpr		+=	en[k].impr;
  average[k][bin].elj		+=	en[k].nbonds - en[k].coulomb;
  average[k][bin].ecoul		+=	en[k].coulomb;
  average[k][bin].epotens	+=  en[k].potens;
  #ifdef SASA
  average[k][bin].esasa		+=  en[k].esasa;
  #endif
  end_to_end(k,SITE1,SITE2);
  average[k][bin].d_nc		+=	ordparam[k].d_nc;
  average[k][bin].hel		+=  ordparam[k].hel;
  average[k][bin].con		+=  ordparam[k].con;
  average[k][bin].con_2		+=  ordparam[k].con_2;
  rad_gyr(k);
  average[k][bin].gyr		+=  ordparam[k].gyr;
  average[k][bin].rmsd		+=  ordparam[k].rmsd;
#endif//DOS

  /* ----------	*/
  /*    MMDOS	*/
  /* ----------	*/
#ifdef MMDOS
  int bin =  (int) ((en[k].totals - sim_dos[k].e_begin)/sim_dos[k].e_width);	
  average[k][bin].count ++;
  average[k][bin].ebond		+=	en[k].bond;
  average[k][bin].ebend		+=	en[k].bend;
#ifndef NEUTRAL
  average[k][bin].eurey		+=	en[k].urey;
#endif
  average[k][bin].etors		+=	en[k].tors;
  average[k][bin].eimpr		+=	en[k].impr;
  average[k][bin].elj		+=	en[k].nbonds - en[k].coulomb;
  average[k][bin].ecoul		+=	en[k].coulomb;
  average[k][bin].epotens	+=  en[k].potens;
  #ifdef SASA
  average[k][bin].esasa		+=  en[k].esasa;
  #endif
  end_to_end(k,SITE1,SITE2);
  average[k][bin].d_nc		+=	ordparam[k].d_nc;
  average[k][bin].hel		+=  ordparam[k].hel;
  average[k][bin].con		+=  ordparam[k].con;
  average[k][bin].con_2		+=  ordparam[k].con_2;
  average[k][bin].rmsd		+=  ordparam[k].rmsd;
#endif//MMDOS

  /* ----------	*/
  /*    CTDOS	*/
  /* ----------	*/
#ifdef CTDOS
  int bin =  (int) ((en[k].potens - sim_dos[k].e_begin)/sim_dos[k].e_width);	
  average[k][bin].count ++;
  average[k][bin].ebond		+=	en[k].bond;
  average[k][bin].ebend		+=	en[k].bend;
#ifndef NEUTRAL
  average[k][bin].eurey		+=	en[k].urey;
#endif
  average[k][bin].etors		+=	en[k].tors;
  average[k][bin].eimpr		+=	en[k].impr;
  average[k][bin].elj		+=	en[k].nbonds - en[k].coulomb;
  average[k][bin].ecoul		+=	en[k].coulomb;
  average[k][bin].epotens	+=  en[k].potens;
  #ifdef SASA
  average[k][bin].esasa		+=  en[k].esasa;
  #endif
  end_to_end(k,SITE1,SITE2);
  average[k][bin].d_nc		+=	ordparam[k].d_nc;
  average[k][bin].hel		+=  ordparam[k].hel;
  average[k][bin].con		+=  ordparam[k].con;
  average[k][bin].con_2		+=  ordparam[k].con_2;
  rad_gyr(k);
  average[k][bin].gyr		+=  ordparam[k].gyr;
  configtemp(k);
  dos_hist[k][bin].ct_num	+= config[k].num;
  dos_hist[k][bin].ct_den	+= config[k].den;
  dos_hist[k][bin].config_T	 = dos_hist[k][bin].ct_num/dos_hist[k][bin].ct_den/RG/0.001;
  average[k][bin].rmsd		+=  ordparam[k].rmsd;
#endif//CTDOS

#ifdef XEDOS
  int bin = (int) ((ordparam[k].d_nc - sim_dos[k].l_begin)/sim_dos[k].l_width);	
  average[k][bin].count ++;
  average[k][bin].ebond		+=	en[k].bond;
  average[k][bin].ebend		+=	en[k].bend;
#ifndef NEUTRAL
  average[k][bin].eurey		+=	en[k].urey;
#endif
  average[k][bin].etors		+=	en[k].tors;
  average[k][bin].eimpr		+=	en[k].impr;
  average[k][bin].elj		+=	en[k].nbonds - en[k].coulomb;
  average[k][bin].ecoul		+=	en[k].coulomb;
  average[k][bin].epotens	+=  en[k].potens;

  average[k][bin].force_1	+=  ordparam[k].force_1;
  average[k][bin].force_2	+=  ordparam[k].force_2;
  
  #ifdef SASA
  average[k][bin].esasa		+=  en[k].esasa;
  #endif
  average[k][bin].d_nc		+=	ordparam[k].d_nc;
  average[k][bin].hel		+=  ordparam[k].hel;
  average[k][bin].con		+=  ordparam[k].con;
  average[k][bin].con_2		+=  ordparam[k].con_2;
  rad_gyr(k);
  average[k][bin].gyr		+=  ordparam[k].gyr;
  average[k][bin].rmsd		+=  ordparam[k].rmsd;
#endif//XEDOS
  /* ----------	*/
  /* FX_EDOS	*/
  /* ----------	*/
#ifdef FX_EDOS
  int bin =  (int) ((en[k].potens - sim_dos[k].e_begin)/sim_dos[k].e_width);
  dos_hist[k][bin].f_of_e = (double)(dos_hist[k][bin].np_of_e) /(double)(dos_hist[k][bin].nw_of_e);
  average[k][bin].count ++;
  average[k][bin].ebond		+=	en[k].bond;
  average[k][bin].ebend		+=	en[k].bend;
#ifndef NEUTRAL
  average[k][bin].eurey		+=	en[k].urey;
#endif
  average[k][bin].etors		+=	en[k].tors;
  average[k][bin].eimpr		+=	en[k].impr;
  average[k][bin].elj		+=	en[k].nbonds - en[k].coulomb;
  average[k][bin].ecoul		+=	en[k].coulomb;
  average[k][bin].epotens	+=  en[k].potens;
  #ifdef SASA
  average[k][bin].esasa		+=  en[k].esasa;
  #endif
  end_to_end(k,SITE1,SITE2);
  average[k][bin].d_nc		+=	ordparam[k].d_nc;
  average[k][bin].hel		+=  ordparam[k].hel;
  average[k][bin].con		+=  ordparam[k].con;
  average[k][bin].con_2		+=  ordparam[k].con_2;
  rad_gyr(k);
  average[k][bin].gyr		+=  ordparam[k].gyr;
  average[k][bin].rmsd		+=  ordparam[k].rmsd;
#endif//FX_EDOS

#ifdef TDXEDOS
  int bin1 = (int) ((ordparam[k].x1 - sim_dos[k].x1_begin)/sim_dos[k].x1_width);
  int bin2 = (int) ((ordparam[k].x2 - sim_dos[k].x2_begin)/sim_dos[k].x2_width);
  average[k][bin1][bin2].count ++;
  average[k][bin1][bin2].ebond		+=	en[k].bond;
  average[k][bin1][bin2].ebend		+=	en[k].bend;
#ifndef NEUTRAL
  average[k][bin1][bin2].eurey		+=	en[k].urey;
#endif
  average[k][bin1][bin2].etors		+=	en[k].tors;
  average[k][bin1][bin2].eimpr		+=	en[k].impr;
  average[k][bin1][bin2].elj		+=	en[k].nbonds - en[k].coulomb;
  average[k][bin1][bin2].ecoul		+=	en[k].coulomb;
  average[k][bin1][bin2].epotens	+=  en[k].potens;

  average[k][bin1][bin2].force_1	+=  ordparam[k].force_1;
  average[k][bin1][bin2].force_2	+=  ordparam[k].force_2;
  
  #ifdef SASA
  average[k][bin1][bin2].esasa		+=  en[k].esasa;
  #endif
  average[k][bin1][bin2].d_nc		+=	ordparam[k].d_nc;
  average[k][bin1][bin2].hel		+=  ordparam[k].hel;
  average[k][bin1][bin2].con		+=  ordparam[k].con;
  average[k][bin].con_2		+=  ordparam[k].con_2;
  rad_gyr(k);
  average[k][bin1][bin2].gyr		+=  ordparam[k].gyr;
  average[k][bin1][bin2].rmsd		+=  ordparam[k].rmsd;
  average[k][bin1][bin2].x1			+=  ordparam[k].x1;
  average[k][bin1][bin2].x2			+=  ordparam[k].x2;
#endif//TDXEDOS

  /* The following added for sampling phi and psi	*/
  for(int i=0; i<torsN[k]; i++) {
	  tors[k][i].phia += tors[k][i].phi;
	  tors[k][i].psia += tors[k][i].psi;
	  tors[k][i].thetaa += tors[k][i].theta;
	  tors[k][i].count ++;
  }
}
