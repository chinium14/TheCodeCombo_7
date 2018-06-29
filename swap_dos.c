#ifdef DOS
/* ======================================================================== */
/* swap_dos.cpp                                                             */
/*                                                                          */
/*		This subroutine swap the conformations in adjacent boxes for density*/
/* pf states simulation. the only criterion for accepting a proposed swap   */
/* is that the two conformations should have energy in the overlapping range*/
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

#ifndef MPI
void nblist		  (int);
double ran2		  (void);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
void flat_histogram(int);
void merge_gofe	  (int);	

int swap_dos(unsigned long iter, int ibox)
{
	int k=ibox;
	double e1 = en[k].potens;
	double e2 = en[k+1].potens;
	int accep_crit=0;
	/*--------------------------------------------------------------------------*/
	/* Determine the old and new bins of each box.                              */
	/*--------------------------------------------------------------------------*/
	int old_bin_1			=  (int) ((e1 - sim_dos[k].e_begin)/sim_dos[k].e_width);
	int new_bin_1			=  (int) ((e2 - sim_dos[k].e_begin)/sim_dos[k].e_width);
	int old_bin_2			=  (int) ((e2 - sim_dos[k+1].e_begin)/sim_dos[k+1].e_width);
	int new_bin_2			=  (int) ((e1 - sim_dos[k+1].e_begin)/sim_dos[k+1].e_width);

	/*--------------------------------------------------------------------------*/
	/* Calculate the accpetance criteria.                                       */
	/*--------------------------------------------------------------------------*/
	/*-------------------------------------------*/
	/* Check to make sure sturctures are in      */
	/* overlap region.  Then check the           */
	/* acceptance criteria.                      */
	/*-------------------------------------------*/
    if ((e2 < sim_dos[k].e_begin || e2 >= sim_dos[k].e_end) || (e1 < sim_dos[k+1].e_begin || e1 >= sim_dos[k+1].e_end))
		accep_crit = 0;
	else{
		if(exp(dos_hist[k][old_bin_1].g_of_e-dos_hist[k][new_bin_1].g_of_e+dos_hist[k+1][old_bin_2].g_of_e-dos_hist[k+1][new_bin_2].g_of_e) > ran2())
			accep_crit = 1;
		else
			accep_crit = 0;
	}

//	if ((e2 >= sim_dos[k].e_begin && e2 < sim_dos[k].e_end) && (e1 >= sim_dos[k+1].e_begin && e1 < sim_dos[k+1].e_end)) 
//		accep_crit= 1;


		if(accep_crit==1)
		{
			for(int i =0; i< box[k].boxns; i++){
				atom_temp[k][i]		= atom[k][i];				/* Back up coordinates of box k			*/
				atnopbc_temp[k][i]	= atnopbc[k][i];			/* Back up coordinates of box k			*/
				atom[k][i]			= atom[k+1][i];				/* Transfer coordinates of k+1 to k		*/
				atnopbc[k][i]		= atnopbc[k+1][i];			/* Transfer coordinates of k+1 to k		*/
				atom[k+1][i]		= atom_temp[k][i];			/* Transfer coordinates of k to k+1		*/
				atnopbc[k+1][i]		= atnopbc_temp[k][i];		/* Transfer coordinates of k to k+1		*/
				ff_temp[k][i]		= ff[k][i];
				ff[k][i]			= ff[k+1][i];
				ff[k+1][i]			= ff_temp[k][i];
				vv_temp[k][i]		= vv[k][i];
				vv[k][i]			= vv[k+1][i];
				vv[k+1][i]			= vv_temp[k][i];
				uu[k][i]			= vv[k][i];
				uu[k+1][i]			= vv[k+1][i];
			}
				
				/*======================================*/
				/* Tranfer data from box k to temp		*/
				/*======================================*/
				en_temp[k]			= en[k];					
				
				
				/*======================================*/
				/* Tranfer data from box k+1 to box k	*/
				/*======================================*/

				en[k]				= en[k+1];									
				#ifdef NLIST
					 nblist(k);									
				#endif
				/*======================================*/
				/* Tranfer data from temp to box k+1	*/
				/*======================================*/

				en[k+1]				= en_temp[k];						
				#ifdef NLIST
					 nblist(k+1);								
				#endif
	/*--------------------------------------------------------------------------*/
	/* The following added to take into account the effect of swapping			*/
	/*--------------------------------------------------------------------------*/
			dos_hist[k][new_bin_1].h_of_e ++;
			dos_hist[k][new_bin_1].g_of_e += sim_dos[k].mod_f;
			flat_histogram(k);
			dos_hist[k+1][new_bin_2].h_of_e ++;
			dos_hist[k+1][new_bin_2].g_of_e += sim_dos[k+1].mod_f;
			flat_histogram(k+1);
	/*--------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------*/
			if (sim_dos[k].mod_f < MERGE_F && sim_dos[k+1].mod_f < MERGE_F && 
				((iter%(MERGE_N*(unsigned long)sim_hyb.cyc_swap)==0) || (iter%(((unsigned long)(MERGE_N+1))*(unsigned long)sim_hyb.cyc_swap)==0))) merge_gofe(k);
			return 1;
		}
		else {
	/*--------------------------------------------------------------------------*/
	/* The following added to take into account the effect of swapping			*/
	/*--------------------------------------------------------------------------*/
			dos_hist[k][old_bin_1].h_of_e ++;
			dos_hist[k][old_bin_1].g_of_e += sim_dos[k].mod_f;
			flat_histogram(k);
			dos_hist[k+1][old_bin_2].h_of_e ++;
			dos_hist[k+1][old_bin_2].g_of_e += sim_dos[k+1].mod_f;
			flat_histogram(k+1);
	/*--------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------*/
			if (sim_dos[k].mod_f < MERGE_F && sim_dos[k+1].mod_f < MERGE_F && 
				((iter%(MERGE_N*(unsigned long)sim_hyb.cyc_swap)==0) || (iter%(((unsigned long)(MERGE_N+1))*(unsigned long)sim_hyb.cyc_swap)==0))) merge_gofe(k);
			return 0;
		}
		
}

void merge_gofe(int k){
	/* Calculate the number of overlapping bins between the two windows	*/

	int o_bins		= (int)((sim_dos[k].e_end - sim_dos[k+1].e_begin)/sim_dos[k].e_width);
/*
	// This version was for simple merging of g_of_e
  double sum=0.0;
	for(int i=0; i<o_bins; i++){
		double diff = dos_hist[k][sim_dos[k].e_bins - o_bins +i].g_of_e - dos_hist[k+1][i].g_of_e;
		sum			+= diff;
	}
	sum /= o_bins;
	for(int j=0; j< sim_dos[k+1].e_bins; j++){
		dos_hist[k+1][j].g_of_e += sum;
	}
	for (int i=0; i<o_bins; i++){
		dos_hist[k+1][i].g_of_e = (dos_hist[k+1][i].g_of_e + dos_hist[k][sim_dos[k].e_bins - o_bins +i].g_of_e)*0.5;
		dos_hist[k][sim_dos[k].e_bins - o_bins +i].g_of_e = dos_hist[k+1][i].g_of_e;
	}
*/

	// This version merges g_of_e based on Temperature information
	
	double T_lwin[NBINS];
	double T_rwin[NBINS];
	double T_mean[NBINS];
	/* j keep track of indices for the left window	*/
	int j = sim_dos[k].e_bins - o_bins;
	
	/* ----------------------------------------	*/
	/* First take the differential of g_of_e in	*/
	/* both windows and average them out		*/
	/* ----------------------------------------	*/
	
	for(int i=0; i<o_bins; i++){
		if (i!=0 && i!=  o_bins-1){
			T_lwin[i] =  dos_hist[k][j+i+1].g_of_e - dos_hist[k][j+i-1].g_of_e;
			T_rwin[i] =  dos_hist[k+1][i+1].g_of_e - dos_hist[k+1][i-1].g_of_e;
		}
		else if (i==0){ // First element of right window
			T_lwin[i] =  dos_hist[k][j+i+1].g_of_e - dos_hist[k][j+i-1].g_of_e;
			T_rwin[i] =	 dos_hist[k+1][2].g_of_e - dos_hist[k+1][0].g_of_e;

		}
		else {			// Last element of left window
			T_lwin[i] =	 T_lwin[i-1];
			T_rwin[i] =  dos_hist[k+1][i+1].g_of_e - dos_hist[k+1][i-1].g_of_e;
		}
		T_mean[i]	  =  (T_lwin[i]+T_rwin[i])/4.0; // 4 takes care of average and twice binwidth
	}

	/* ----------------------------------------	*/
	/* Now first adjust the g_of_e for the left	*/
	/* window in the overlapping range. 		*/
	/* ----------------------------------------	*/

	for(int i=0; i<o_bins; i++){
		double g_old = dos_hist[k][j+i].g_of_e; //added for H(E)
		dos_hist[k][j+i].g_of_e =  dos_hist[k][j+i-1].g_of_e + T_mean[i];
		dos_hist[k][j+i].h_of_e += (int)((dos_hist[k][j+i].g_of_e - g_old)/sim_dos[k].mod_f); //added for H(E)
		if(dos_hist[k][j+i].h_of_e <0) dos_hist[k][j+i].h_of_e = 0;
	}
	/* ----------------------------------------	*/
	/* Now adjust the g_of_e for the right win 	*/
	/* outside the overlapping range:just shift	*/
	/* ----------------------------------------	*/
	double shift = dos_hist[k+1][o_bins-1].g_of_e -dos_hist[k][j+o_bins-1].g_of_e;
	for(int i=o_bins; i< sim_dos[k+1].e_bins; i++){
		dos_hist[k+1][i].g_of_e -= shift;
	}
	/* ----------------------------------------	*/
	/* Now adjust the g_of_e for the right win	*/
	/* in the overlapping range.		 		*/
	/* ----------------------------------------	*/
	for(int i=0; i<o_bins; i++){
		double g_old = dos_hist[k+1][i].g_of_e; // added for H(E)
		dos_hist[k+1][i].g_of_e =  dos_hist[k][j+i].g_of_e;	
		dos_hist[k+1][i].h_of_e +=  (int)((dos_hist[k+1][i].g_of_e - g_old + shift)/sim_dos[k+1].mod_f);
		if(dos_hist[k+1][i].h_of_e <0) dos_hist[k+1][i].h_of_e = 0;
	}
}
#endif
#endif
