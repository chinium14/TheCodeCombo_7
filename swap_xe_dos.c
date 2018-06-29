/* ======================================================================== */
/* swap_xe_dos.cpp                                                          */
/*                                                                          */
/*		This subroutine swap the conformations in adjacent boxes for density*/
/* pf states simulation. the only criterion for accepting a proposed swap   */
/* is that the two conformations should have energy in the overlapping range*/
/*                                                                          */
/* ======================================================================== */

#include "defines.h"
#ifdef XEDOS
#ifndef MPI
void nblist		  (int);
double ran2		  (void);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
void flat_histogram(int);
void merge_gofl	  (int);
void end_to_end	  (int,int,int);	

int swap_xe_dos(unsigned long iter, int ibox)
{
	int k=ibox;
	double l1 = ordparam[k].d_nc;
	double l2 = ordparam[k+1].d_nc;
	int accep_crit=0;
	/*--------------------------------------------------------------------------*/
	/* The following added to take into account the effect of swapping			*/
	/*--------------------------------------------------------------------------*/
	int old_bin_1			=  (int) ((l1 - sim_dos[k].l_begin)/sim_dos[k].l_width);
	int new_bin_1			=  (int) ((l2 - sim_dos[k].l_begin)/sim_dos[k].l_width);
	int old_bin_2			=  (int) ((l2 - sim_dos[k+1].l_begin)/sim_dos[k+1].l_width);
	int new_bin_2			=  (int) ((l1 - sim_dos[k+1].l_begin)/sim_dos[k+1].l_width);
	/*--------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------*/
	if ((l2 < sim_dos[k].l_begin || l2 >= sim_dos[k].l_end) || (l1  < sim_dos[k+1].l_begin || l1 >= sim_dos[k+1].l_end))
		accep_crit = 0;
	else{
		if(exp(dos_hist[k][old_bin_1].g_of_l-dos_hist[k][new_bin_1].g_of_l+dos_hist[k+1][old_bin_2].g_of_l-dos_hist[k+1][new_bin_2].g_of_l) > ran2())
			accep_crit = 1;
		else
			accep_crit = 0;
	}



//	if ((l2 >= sim_dos[k].l_begin && l2 < sim_dos[k].l_end) && (l1 >= sim_dos[k+1].l_begin && l1 < sim_dos[k+1].l_end)) 
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
			dos_hist[k][new_bin_1].h_of_l ++;
			dos_hist[k][new_bin_1].g_of_l += sim_dos[k].mod_f;
			flat_histogram(k);
			dos_hist[k+1][new_bin_2].h_of_l ++;
			dos_hist[k+1][new_bin_2].g_of_l += sim_dos[k+1].mod_f;
			flat_histogram(k+1);
	/*--------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------*/
			if (sim_dos[k].mod_f < MERGE_F && sim_dos[k+1].mod_f < MERGE_F && 
				((iter%(MERGE_N*(unsigned long)sim_hyb.cyc_swap)==0) || (iter%(((unsigned long)(MERGE_N+1))*(unsigned long)sim_hyb.cyc_swap)==0))) merge_gofl(k);
			
			end_to_end(k,SITE1,SITE2);
			end_to_end(k+1,SITE1,SITE2);
			return 1;
		}
		else {
	/*--------------------------------------------------------------------------*/
	/* The following added to take into account the effect of swapping			*/
	/*--------------------------------------------------------------------------*/
			dos_hist[k][old_bin_1].h_of_l ++;
			dos_hist[k][old_bin_1].g_of_l += sim_dos[k].mod_f;
			flat_histogram(k);
			dos_hist[k+1][old_bin_2].h_of_l ++;
			dos_hist[k+1][old_bin_2].g_of_l += sim_dos[k+1].mod_f;
			flat_histogram(k+1);
	/*--------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------*/
			if (sim_dos[k].mod_f < MERGE_F && sim_dos[k+1].mod_f < MERGE_F && 
				((iter%(MERGE_N*(unsigned long)sim_hyb.cyc_swap)==0) || (iter%(((unsigned long)(MERGE_N+1))*(unsigned long)sim_hyb.cyc_swap)==0))) merge_gofl(k);
			
			end_to_end(k,SITE1,SITE2);
			end_to_end(k+1,SITE1,SITE2);
			return 0;
		}
		
}

void merge_gofl(int k){
	/* Calculate the number of overlapping bins between the two windows	*/

	int o_bins		= (int)((sim_dos[k].l_end - sim_dos[k+1].l_begin)/sim_dos[k].l_width);

	// This version is for simple merging of g_of_l
  double sum=0.0;
	for(int i=0; i<o_bins; i++){
		double diff = dos_hist[k][sim_dos[k].l_bins - o_bins +i].g_of_l - dos_hist[k+1][i].g_of_l;
		sum			+= diff;
	}
	sum /= o_bins;
	for(int j=0; j< sim_dos[k+1].l_bins; j++){
		dos_hist[k+1][j].g_of_l += sum;
	}
	for (int i=0; i<o_bins; i++){
		dos_hist[k+1][i].g_of_l = (dos_hist[k+1][i].g_of_l + dos_hist[k][sim_dos[k].l_bins - o_bins +i].g_of_l)*0.5;
		dos_hist[k][sim_dos[k].l_bins - o_bins +i].g_of_l = dos_hist[k+1][i].g_of_l;
	}

}
#endif
#endif
