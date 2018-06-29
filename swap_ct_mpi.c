#ifdef MPI
#ifdef CTDOS
/* ========================================================================= */
/* swap_ct_mpi.cpp                                                              */
/*                                                                           */
/*		This subroutine swap the conformations in adjacent boxes for density */
/* of states simulation. the only criterion for accepting a proposed swap    */
/* is that the two conformations should have energy in the overlapping range */
/* Differnt processes are executed on different processors.                  */
/* ========================================================================= */

#include "defines.h"

void nblist		   (int);
double ran2		   (void);
void flat_histogram(int);
void merge_gofe	   (int, double, double, double);
#ifdef STATUS
void curr_status   (int,int);
#endif

int swap_ct_mpi(unsigned long iter, int flag, int ibox)
{
	int k =ibox;
	/*------------------------------------------------------*/
	/* The energy of this process is saved in e1 whereas	*/
	/* The energy of the partner process is saved in e2		*/
	/*------------------------------------------------------*/
	double e1 = en[k].potens;
	int old_bin_1			=  (int) ((e1 - sim_dos[k].e_begin)/sim_dos[k].e_width);
	MPI_Status status;
#ifdef STATUS
			curr_status(k,5);
#endif

	if(mpi.my_rank!=mpi.p-1 && flag ==1)
	{

		struct msg_swap_request swap_request;
		swap_request.potens		=	en[k].potens;
		swap_request.e_begin	=	sim_dos[k].e_begin;
		swap_request.e_end		=	sim_dos[k].e_end;
		swap_request.mod_f		=	sim_dos[k].mod_f;

		struct msg_swap_accept swap_accept;
		
		/* ------------------------------------------------ */
		/* If I am not the last process and if my flag ==1  */
		/* then send data/swap request to my_rank +1		*/
		/* ------------------------------------------------ */
		MPI_Send(&swap_request,sizeof(swap_request),MPI_CHAR,mpi.my_rank+1,1,MPI_COMM_WORLD);

		/* ------------------------------------------------ */
		/* Receive the information whether move is accepted */
		/* or rejected.										*/
		/* ------------------------------------------------ */

		MPI_Recv(&swap_accept, sizeof(swap_accept), MPI_CHAR,mpi.my_rank+1,2,MPI_COMM_WORLD,&status);

		
		/* ------------------------------------------------ */
		/* If move is rejected no more communication needed */
		/* ------------------------------------------------ */
		if(swap_accept.accept==0){
			dos_hist[k][old_bin_1].h_of_e ++;
			dos_hist[k][old_bin_1].g_of_e += sim_dos[k].mod_f;
			flat_histogram(k); //fprintf(stdout,"box %d and box %d swap rejected \n", mpi.my_rank, mpi.my_rank+1);	 
//			if(sim_dos[k].mod_f <MERGE_F && ((iter%(MERGE_N*sim_hyb.cyc_swap)==0) || (iter%(((int)(MERGE_N*0.5+1))*sim_hyb.cyc_swap)==0))){
			if(sim_dos[k].mod_f <MERGE_F && ((iter%(MERGE_N*(unsigned long)sim_hyb.cyc_swap)==0) || (iter%(((unsigned long)(MERGE_N+1))*(unsigned long)sim_hyb.cyc_swap)==0))){
				MPI_Send(dos_hist[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank+1,14,MPI_COMM_WORLD);
				MPI_Recv(dos_hist[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank+1,15,MPI_COMM_WORLD,&status);
			}

			return 0;
		}// is swap not accepted then for the sending box

		/* ------------------------------------------------ */
		/* If move is accepted then create a back up of		*/
		/* of current structures and then receive the data	*/
		/* from the partner process and after that send the */
		/* back up info to the partner to complete the swap */
		/* ------------------------------------------------ */
		else if (swap_accept.accept==1){
			for(int i =0; i< box[k].boxns; i++){
				atom_temp[k][i]		= atom[k][i];				/* Back up coordinates of box k			*/
				atnopbc_temp[k][i]	= atnopbc[k][i];			/* Back up coordinates of box k			*/
				ff_temp[k][i]		= ff[k][i];
				vv_temp[k][i]		= vv[k][i];
			}
							
				/*======================================*/
				/* Tranfer data from box k to temp		*/
				/*======================================*/
				en_temp[k]			= en[k];


			MPI_Recv(atom[k],   box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank+1,3,MPI_COMM_WORLD,&status);
			MPI_Recv(atnopbc[k],box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank+1,4,MPI_COMM_WORLD,&status);
			MPI_Recv(ff[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank+1,5,MPI_COMM_WORLD,&status);
			MPI_Recv(vv[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank+1,6,MPI_COMM_WORLD,&status);
			MPI_Recv(&en[k],     sizeof(struct energy),MPI_CHAR,mpi.my_rank+1,11,MPI_COMM_WORLD,&status);

			uu[k]			= vv[k];
			en[k].potens	= swap_accept.potens;
			int new_bin_1	= (int) ((en[k].potens - sim_dos[k].e_begin)/sim_dos[k].e_width);

			MPI_Send(atom_temp[k],   box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank+1, 7,MPI_COMM_WORLD);
			MPI_Send(atnopbc_temp[k],box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank+1, 8,MPI_COMM_WORLD);
			MPI_Send(ff_temp[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank+1, 9,MPI_COMM_WORLD);
			MPI_Send(vv_temp[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank+1,10,MPI_COMM_WORLD);
			MPI_Send(&en_temp[k],     sizeof(struct energy),MPI_CHAR,mpi.my_rank+1,12,MPI_COMM_WORLD);
			
		

			dos_hist[k][new_bin_1].h_of_e ++;
			dos_hist[k][new_bin_1].g_of_e += sim_dos[k].mod_f;
			flat_histogram(k);
			#ifdef NLIST
				 nblist(k);									
			#endif
//			if(sim_dos[k].mod_f <MERGE_F && ((iter%(MERGE_N*sim_hyb.cyc_swap)==0) || (iter%(((int)(MERGE_N*0.5+1))*sim_hyb.cyc_swap)==0))){
				 
			if(sim_dos[k].mod_f <MERGE_F && ((iter%(MERGE_N*(unsigned long)sim_hyb.cyc_swap)==0) || (iter%(((unsigned long)(MERGE_N+1))*(unsigned long)sim_hyb.cyc_swap)==0))){
				 MPI_Send(dos_hist[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank+1,14,MPI_COMM_WORLD);
				 MPI_Recv(dos_hist[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank+1,15,MPI_COMM_WORLD,&status);
			}
			return 1;
		}// if swap accepted then for the sending box
	 }// end of mpi.rank !=p-1
	 
	 else if(mpi.my_rank!=0 && flag==0)
	 {
			/* Recv data from my_rank-1 */
		struct msg_swap_request swap_request;
//		swap_request.potens		=	en[k].potens;
//		swap_request.e_begin	=	sim_dos[k].e_begin;
//		swap_request.e_end		=	sim_dos[k].e_end;

		struct msg_swap_accept swap_accept;
		/* Recv data from my_rank-1 */
		MPI_Recv(&swap_request,sizeof(swap_request),MPI_CHAR,mpi.my_rank-1,1,MPI_COMM_WORLD,&status);
		double e2 = swap_request.potens;
		double mod_f = swap_request.mod_f;

		/* ----------------------------------------	*/
		/*	This onwards added later for on the fly	*/
		/*	merging of g_of_e						*/
		double e_begin	=	swap_request.e_begin;
		double e_end	=	swap_request.e_end;
		/* ----------------------------------------	*/

		if ((e2 >= sim_dos[k].e_begin && e2 < sim_dos[k].e_end) && (e1 >= swap_request.e_begin && e1 < swap_request.e_end)) 
		{
			swap_accept.accept = 1;
			swap_accept.potens = e1;
			MPI_Send(&swap_accept,sizeof(swap_accept),   MPI_CHAR,mpi.my_rank-1,2,MPI_COMM_WORLD);
			MPI_Send(atom[k],   box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank-1,3,MPI_COMM_WORLD);
			MPI_Send(atnopbc[k],box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank-1,4,MPI_COMM_WORLD);
			MPI_Send(ff[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank-1,5,MPI_COMM_WORLD);
			MPI_Send(vv[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank-1,6,MPI_COMM_WORLD);
			MPI_Send(&en[k],     sizeof(struct energy),MPI_CHAR,mpi.my_rank-1,11,MPI_COMM_WORLD);
			en[k].potens	= e2;
			int new_bin_1	= (int) ((en[k].potens - sim_dos[k].e_begin)/sim_dos[k].e_width);
			MPI_Recv(atom[k],   box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank-1, 7,MPI_COMM_WORLD,&status);
			MPI_Recv(atnopbc[k],box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank-1, 8,MPI_COMM_WORLD,&status);
			MPI_Recv(ff[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank-1, 9,MPI_COMM_WORLD,&status);
			MPI_Recv(vv[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank-1,10,MPI_COMM_WORLD,&status);
			MPI_Recv(&en[k],     sizeof(struct energy),MPI_CHAR,mpi.my_rank-1,12,MPI_COMM_WORLD,&status);
			uu[k]	=	vv[k];
			dos_hist[k][new_bin_1].h_of_e ++;
			dos_hist[k][new_bin_1].g_of_e += sim_dos[k].mod_f;
			flat_histogram(k);
			#ifdef NLIST
					 nblist(k);									
			#endif
//			if(mod_f<MERGE_F && ((iter%(MERGE_N*sim_hyb.cyc_swap)==0) || (iter%(((int)(MERGE_N*0.5+1))*sim_hyb.cyc_swap)==0))){
			
			if(mod_f<MERGE_F && ((iter%(MERGE_N*(unsigned long)sim_hyb.cyc_swap)==0) || (iter%(((unsigned long)(MERGE_N+1))*(unsigned long)sim_hyb.cyc_swap)==0))){
				MPI_Recv(dos_hist_temp[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank-1,14,MPI_COMM_WORLD,&status);
				if (sim_dos[k].mod_f <MERGE_F) merge_gofe(k,e_begin,e_end, mod_f);
				MPI_Send(dos_hist_temp[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank-1,15,MPI_COMM_WORLD);
			}
			return 1;
		}//end of accept for receiving box
		else {
			swap_accept.accept=0;
			swap_accept.potens=e2;// e2 belongs to myrank-1
			MPI_Send(&swap_accept,sizeof(swap_accept),   MPI_CHAR,mpi.my_rank-1,2,MPI_COMM_WORLD);

			dos_hist[k][old_bin_1].h_of_e ++;
			dos_hist[k][old_bin_1].g_of_e += sim_dos[k].mod_f;
			flat_histogram(k);
//			if(mod_f<MERGE_F && ((iter%(MERGE_N*sim_hyb.cyc_swap)==0) || (iter%(((int)(MERGE_N*0.5+1))*sim_hyb.cyc_swap)==0))){
			if(mod_f<MERGE_F && ((iter%(MERGE_N*(unsigned long)sim_hyb.cyc_swap)==0) || (iter%(((unsigned long)(MERGE_N+1))*(unsigned long)sim_hyb.cyc_swap)==0))){
				MPI_Recv(dos_hist_temp[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank-1,14,MPI_COMM_WORLD,&status);
				if (sim_dos[k].mod_f <MERGE_F) merge_gofe(k,e_begin,e_end,mod_f);
				MPI_Send(dos_hist_temp[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank-1,15,MPI_COMM_WORLD);
			}
			
			return 0;
		}// end of reject for receiving box
	 }// end of rank!=0 loop
	 return 2;
}

void merge_gofe(int k, double e_begin, double e_end, double mod_f){
	int e_bins		= (int)((e_end - e_begin)/sim_dos[k].e_width+0.5);// Number of bins in my_rank -1 histogram
	int o_bins		= (int)((e_end - sim_dos[k].e_begin)/sim_dos[k].e_width);
/*
	double sum=0.0;
	for(int i=0; i<o_bins; i++){
		double diff = dos_hist_temp[k][e_bins - o_bins +i].g_of_e - dos_hist[k][i].g_of_e;
		sum			+= diff;
	}
	sum /= o_bins;
	for(int j=0; j< sim_dos[k].e_bins; j++){
		dos_hist[k][j].g_of_e += sum;
	}
	for (int i=0; i<o_bins; i++){
		dos_hist[k][i].g_of_e = (dos_hist[k][i].g_of_e + dos_hist_temp[k][e_bins - o_bins +i].g_of_e)*0.5;
		dos_hist_temp[k][e_bins - o_bins +i].g_of_e = dos_hist[k][i].g_of_e;
	}
*/
		// This version merges g_of_e based on Temperature information
	
	double T_lwin[NBINS];
	double T_rwin[NBINS];
	double T_mean[NBINS];
	/* j keep track of indices for the left window	*/
	int j = e_bins - o_bins;
	
	/* ----------------------------------------	*/
	/* First take the differential of g_of_e in	*/
	/* both windows and average them out		*/
	/* ----------------------------------------	*/
	
	for(int i=0; i<o_bins; i++){
		if (i!=0 && i!=  o_bins-1){
			T_lwin[i] =  dos_hist_temp[k][j+i+1].g_of_e - dos_hist_temp[k][j+i-1].g_of_e;
			T_rwin[i] =  dos_hist[k][i+1].g_of_e		- dos_hist[k][i-1].g_of_e;
		}
		else if (i==0){ // First element of right window
			T_lwin[i] =  dos_hist_temp[k][j+i+1].g_of_e - dos_hist_temp[k][j+i-1].g_of_e;
			T_rwin[i] =	 dos_hist[k][2].g_of_e	- dos_hist[k][0].g_of_e;

		}
		else {			// Last element of left window
			T_lwin[i] =	 T_lwin[i-1];
			T_rwin[i] =  dos_hist[k][i+1].g_of_e - dos_hist[k][i-1].g_of_e;
		}
		T_mean[i]	  =  (T_lwin[i]+T_rwin[i])/4.0; // 4 takes care of average and twice binwidth
	}

	/* ----------------------------------------	*/
	/* Now first adjust the g_of_e for the left	*/
	/* window in the overlapping range. 		*/
	/* ----------------------------------------	*/

	for(int i=0; i<o_bins; i++){
		double g_old = dos_hist_temp[k][j+i].g_of_e; //added for H(E)
		dos_hist_temp[k][j+i].g_of_e =  dos_hist_temp[k][j+i-1].g_of_e + T_mean[i];
		dos_hist_temp[k][j+i].h_of_e += (int)((dos_hist_temp[k][j+i].g_of_e - g_old)/mod_f); //added for H(E)
		if(dos_hist_temp[k][j+i].h_of_e <0) dos_hist_temp[k][j+i].h_of_e = 0;
	}
	/* ----------------------------------------	*/
	/* Now adjust the g_of_e for the right win 	*/
	/* outside the overlapping range:just shift	*/
	/* ----------------------------------------	*/
	double shift =	dos_hist[k][o_bins-1].g_of_e - dos_hist_temp[k][j+o_bins-1].g_of_e;
	for(int i=o_bins; i< sim_dos[k].e_bins; i++){
		dos_hist[k][i].g_of_e -= shift;
	}
	/* ----------------------------------------	*/
	/* Now adjust the g_of_e for the right win	*/
	/* in the overlapping range.		 		*/
	/* ----------------------------------------	*/
	for(int i=0; i<o_bins; i++){
		double g_old = dos_hist[k][i].g_of_e; // added for H(E)
		dos_hist[k][i].g_of_e =  dos_hist_temp[k][j+i].g_of_e;
		dos_hist[k][i].h_of_e +=  (int)((dos_hist[k][i].g_of_e - g_old + shift)/sim_dos[k].mod_f);
		if(dos_hist[k][i].h_of_e <0) dos_hist[k][i].h_of_e = 0;
	}
		/* ----------------------------------------	*/
	/* Now merge the numberator and denominator */
	/* of config T in the overlapping range.	*/
	/* ----------------------------------------	*/
	for(int i=0; i<o_bins; i++){
		dos_hist_temp[k][j+i].ct_num =  dos_hist_temp[k][j+i].ct_num + dos_hist[k][i].ct_num;
		dos_hist_temp[k][j+i].ct_den =  dos_hist_temp[k][j+i].ct_den + dos_hist[k][i].ct_den;
		dos_hist[k][i].ct_num =  dos_hist_temp[k][j+i].ct_num;
		dos_hist[k][i].ct_den =  dos_hist_temp[k][j+i].ct_den;
	}
}
#endif//DOS
#endif//MPI
