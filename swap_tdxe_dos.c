/* ======================================================================== */
/* swap_tdxe_dos.cpp                                                          */
/*                                                                          */
/*		This subroutine swap the conformations in adjacent boxes for density*/
/* pf states simulation. the only criterion for accepting a proposed swap   */
/* is that the two conformations should have energy in the overlapping range*/
/*                                                                          */
/* ======================================================================== */

#include "defines.h"
#ifdef TDXEDOS
#ifndef MPI
void nblist		  (int);
double ran2		  (void);
int ran_int(int, int);
void flat_histogram_td(int);
void merge_gofl	  (int);
void end_to_end	  (int,int,int);	
double wall_angle(int, int);
double dos_interp_2 (int, double, double);
int swap_tdxe_dos(int box_1,int box_2)
{
	int box1 = box_1;
	int box2 = box_2;
    /*------------------------------------------*/
    /*Determine the values of the reaction      */
    /*coordiates for each box.                  */
    /*------------------------------------------*/
	end_to_end(box1,SITE1,SITE2);
	end_to_end(box2,SITE1,SITE2);
	double x1_box1	=	ordparam[box1].d_nc;
	double x1_box2	=	ordparam[box2].d_nc;
	double x2_box1	=	wall_angle(box1,wall[box1].angle_site);
	double x2_box2	=	wall_angle(box2,wall[box2].angle_site);
	
	int accep_crit=0;
	/*--------------------------------------------------------------------------*/
	/* The following added to take into account the effect of swapping			*/
	/*--------------------------------------------------------------------------*/
	int old_bin_x1_box1	=  (int) ((x1_box1 - sim_dos[box1].x1_begin)/sim_dos[box1].x1_width);
	int old_bin_x2_box1 =  (int) ((x2_box1 - sim_dos[box1].x2_begin)/sim_dos[box1].x2_width);
	int new_bin_x1_box1	=  (int) ((x1_box2 - sim_dos[box1].x1_begin)/sim_dos[box1].x1_width);
	int new_bin_x2_box1 =  (int) ((x2_box2 - sim_dos[box1].x2_begin)/sim_dos[box1].x2_width);

	int old_bin_x1_box2 =  (int) ((x1_box2 - sim_dos[box2].x1_begin)/sim_dos[box2].x1_width);
	int old_bin_x2_box2	=  (int) ((x2_box2 - sim_dos[box2].x2_begin)/sim_dos[box2].x2_width);
	int new_bin_x1_box2 =  (int) ((x1_box1 - sim_dos[box2].x1_begin)/sim_dos[box2].x1_width);
    int new_bin_x2_box2 =  (int) ((x2_box1 - sim_dos[box2].x2_begin)/sim_dos[box2].x2_width);

    /*--------------------------------------------------------------------------*/
    /* To see if the proposed swap is accepted, first see if the structures lie */
    /* in the overlapping range, and then evaluate the acceptance criterion.    */
    /*--------------------------------------------------------------------------*/
	if ((x1_box2 >= sim_dos[box1].x1_begin && x1_box2 < sim_dos[box1].x1_end) && 
		(x2_box2 >= sim_dos[box1].x2_begin && x2_box2 < sim_dos[box1].x2_end) &&
		(x1_box1 >= sim_dos[box2].x1_begin && x1_box1 < sim_dos[box2].x1_end) &&
		(x2_box1 >= sim_dos[box2].x2_begin && x2_box1 < sim_dos[box2].x2_end)){
        /* ------------------------------------------------ */
        /* Determine the new and old density of states by   */
        /* 2D interpolation.                                */
        /* ------------------------------------------------ */
        double g_old_box1, g_new_box1, g_old_box2, g_new_box2;

        g_old_box1 = dos_interp_2(box1,x1_box1,x2_box1);
        g_new_box1 = dos_interp_2(box1,x1_box2,x2_box2);
		g_old_box2 = dos_interp_2(box2,x1_box2,x2_box2);
		g_new_box2 = dos_interp_2(box2,x1_box1,x2_box1);
    	
		double arg = g_old_box1 - g_new_box1 + g_old_box2 - g_new_box2;
		if(arg>0.0) accep_crit = 1;
		else if(exp(arg) > ran2()) accep_crit = 1;
		else accep_crit = 0;
	}else//rejected because out of bounds 
		accep_crit = 0;

	/*======================================*/
	/* Transfer the data between each box.  */
	/*======================================*/
	if(accep_crit==1){
		for(int i =0; i< box[box1].boxns; i++){
			atom_temp[box1][i]		= atom[box1][i];				/* Back up coordinates of box k			*/
			atnopbc_temp[box1][i]	= atnopbc[box1][i];				/* Back up coordinates of box k			*/
			atom[box1][i]			= atom[box2][i];				/* Transfer coordinates of k+1 to k		*/
			atnopbc[box1][i]		= atnopbc[box2][i];				/* Transfer coordinates of k+1 to k		*/
			atom[box2][i]			= atom_temp[box1][i];			/* Transfer coordinates of k to k+1		*/
			atnopbc[box2][i]		= atnopbc_temp[box1][i];		/* Transfer coordinates of k to k+1		*/
			ff_temp[box1][i]		= ff[box1][i];
			ff[box1][i]				= ff[box2][i];
			ff[box2][i]				= ff_temp[box1][i];
			vv_temp[box1][i]		= vv[box1][i];
			vv[box1][i]				= vv[box2][i];
			vv[box2][i]				= vv_temp[box1][i];
			uu[box1][i]				= vv[box1][i];
			uu[box2][i]				= vv[box2][i];
		}
				
		en_temp[box1]			= en[box1];					
		en[box1]				= en[box2];									
		en[box2]				= en_temp[box1];						
		#ifdef NLIST
		nblist(box1);									
		nblist(box2);	
		#endif

		/*-----------------------------------------------------------------------*/
		/* Update the histograms and dos estimate.                  		     */
		/*-----------------------------------------------------------------------*/
		dos_hist[box1][new_bin_x1_box1][new_bin_x2_box1].h_of_l ++;
		dos_hist[box1][new_bin_x1_box1][new_bin_x2_box1].g_of_l += sim_dos[box1].mod_f;
		flat_histogram_td(box1);
		dos_hist[box2][new_bin_x1_box2][new_bin_x2_box2].h_of_l ++;
		dos_hist[box2][new_bin_x1_box2][new_bin_x2_box2].g_of_l += sim_dos[box2].mod_f;
		flat_histogram_td(box2);

        /*-----------------------------------------------------------------------*/
        /* Merge the estimates of the density of states.                         */
        /*-----------------------------------------------------------------------*/
//		if (sim_dos[box1].mod_f < MERGE_F && sim_dos[box2].mod_f < MERGE_F && 
//				((iter%(MERGE_N*sim_hyb.cyc_swap)==0) || (iter%(((int)(MERGE_N+1))*sim_hyb.cyc_swap)==0))) merge_gofl(k);

        /*-----------------------------------------------------------------------*/
        /* Update the reaction coordinates.                                      */
        /*-----------------------------------------------------------------------*/
		end_to_end(box1,SITE1,SITE2);
		end_to_end(box2,SITE1,SITE2);
		ordparam[box1].x1 = ordparam[box1].d_nc;
		ordparam[box1].x2 = wall_angle(box1,wall[box1].angle_site);
		ordparam[box2].x1 = ordparam[box2].d_nc;
		ordparam[box2].x2 = wall_angle(box2,wall[box2].angle_site);
		return 1;
	}
	else{
        /*-----------------------------------------------------------------------*/
        /* Update the histograms and dos estimate.                               */
        /*-----------------------------------------------------------------------*/
		dos_hist[box1][old_bin_x1_box1][old_bin_x2_box1].h_of_l ++;
		dos_hist[box1][old_bin_x1_box1][old_bin_x2_box1].g_of_l += sim_dos[box1].mod_f;
		flat_histogram_td(box1);
		dos_hist[box2][old_bin_x1_box2][old_bin_x2_box2].h_of_l ++;
		dos_hist[box2][old_bin_x1_box2][old_bin_x2_box2].g_of_l += sim_dos[box2].mod_f;
		flat_histogram_td(box2);

       	/*-----------------------------------------------------------------------*/
       	/* Merge the estimates of the density of states.                         */
       	/*-----------------------------------------------------------------------*/
//		if (sim_dos[box1].mod_f < MERGE_F && sim_dos[box2].mod_f < MERGE_F && 
//			((iter%(MERGE_N*sim_hyb.cyc_swap)==0) || (iter%(((int)(MERGE_N+1))*sim_hyb.cyc_swap)==0))) merge_gofl(k);

        /*-----------------------------------------------------------------------*/
        /* Update the reaction coordinates.                                      */
        /*-----------------------------------------------------------------------*/
        end_to_end(box1,SITE1,SITE2);
        end_to_end(box2,SITE1,SITE2);
        ordparam[box1].x1 = ordparam[box1].d_nc;
        ordparam[box1].x2 = wall_angle(box1,wall[box1].angle_site);
        ordparam[box2].x1 = ordparam[box2].d_nc;
        ordparam[box2].x2 = wall_angle(box2,wall[box2].angle_site);			
		return 0;
	}
		
}

/*void merge_gofl(int k){

	int o_bins		= (int)((sim_dos[box1].l_end - sim_dos[box2].l_begin)/sim_dos[box1].l_width);

	// This version is for simple merging of g_of_l
  double sum=0.0;
	for(int i=0; i<o_bins; i++){
		double diff = dos_hist[box1][sim_dos[box1].l_bins - o_bins +i].g_of_l - dos_hist[box2][i].g_of_l;
		sum			+= diff;
	}
	sum /= o_bins;
	for(int j=0; j< sim_dos[box2].l_bins; j++){
		dos_hist[box2][j].g_of_l += sum;
	}
	for (int i=0; i<o_bins; i++){
		dos_hist[box2][i].g_of_l = (dos_hist[box2][i].g_of_l + dos_hist[box1][sim_dos[box1].l_bins - o_bins +i].g_of_l)*0.5;
		dos_hist[box1][sim_dos[box1].l_bins - o_bins +i].g_of_l = dos_hist[box2][i].g_of_l;
	}

}*/
#endif
#endif
