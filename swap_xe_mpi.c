/* ========================================================================= */
/* swap_mpi.cpp                                                              */
/*                                                                           */
/*	This subroutine swap the conformations in adjacent boxes for density */
/* of states simulation. the only criterion for accepting a proposed swap    */
/* is that the two conformations should have energy in the overlapping range */
/* Differnt processes are executed on different processors.                  */
/* ========================================================================= */

#include "defines.h"
#ifdef MPI
#ifdef XEDOS
void nblist		  (int);
double ran2		  (void);
void flat_histogram	  (int);
void merge_gofl	   	  (int, double, double, double);
void end_to_end	  	  (int,int,int);


int swap_xe_mpi(unsigned long iter, int flag, int ibox)
{
	int k =ibox;
	/*------------------------------------------------------*/
	/* The energy of this process is saved in e1 whereas	*/
	/* The energy of the partner process is saved in e2	*/
	/*------------------------------------------------------*/
	double l1 = ordparam[k].d_nc;
	
	int old_bin_1			= (int) ((l1 - sim_dos[k].l_begin)/sim_dos[k].l_width);
	MPI_Status status;

	if(mpi.my_rank!=mpi.p-1 && flag ==1)
	{

		struct msg_swap_request swap_request;
		swap_request.potens		=	en[k].potens;
		swap_request.d_nc		=	ordparam[k].d_nc;  
		swap_request.l_begin		=	sim_dos[k].l_begin;
		swap_request.l_end		=	sim_dos[k].l_end;
		swap_request.mod_f		=	sim_dos[k].mod_f;

		struct msg_swap_accept swap_accept;
		
		/* ------------------------------------------------ */
		/* If I am not the last process and if my flag ==1  */
		/* then send data/swap request to my_rank +1	    */
		/* ------------------------------------------------ */
		MPI_Send(&swap_request,sizeof(swap_request),MPI_CHAR,mpi.my_rank+1,1,MPI_COMM_WORLD);

		/* ------------------------------------------------ */
		/* Receive the information whether move is accepted */
		/* or rejected.					    */
		/* ------------------------------------------------ */

		MPI_Recv(&swap_accept, sizeof(swap_accept), MPI_CHAR,mpi.my_rank+1,2,MPI_COMM_WORLD,&status);

		
		/* ------------------------------------------------ */
		/* If move is rejected no more communication needed */
		/* ------------------------------------------------ */
		if(swap_accept.accept==0){
			dos_hist[k][old_bin_1].h_of_l ++;
			dos_hist[k][old_bin_1].g_of_l += sim_dos[k].mod_f;
			flat_histogram(k); 
			if(sim_dos[k].mod_f <MERGE_F && ((iter%(MERGE_N*(unsigned long)sim_hyb.cyc_swap)==0) || 
						(iter%(((unsigned long)(MERGE_N+1))*(unsigned long)sim_hyb.cyc_swap)==0))){
				MPI_Send(dos_hist[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank+1,14,MPI_COMM_WORLD);
				MPI_Recv(dos_hist[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank+1,15,MPI_COMM_WORLD,&status);
			}

			return 0;
		}// is swap not accepted then for the sending box

		/* ------------------------------------------------ */
		/* If move is accepted then create a back up of	    */
		/* of current structures and then receive the data  */
		/* from the partner process and after that send the */
		/* back up info to the partner to complete the swap */
		/* ------------------------------------------------ */
		else if (swap_accept.accept==1){
			for(int i =0; i< box[k].boxns; i++){
				atom_temp[k][i]		= atom[k][i];			/* Back up coordinates of box k	*/
				atnopbc_temp[k][i]	= atnopbc[k][i];		/* Back up coordinates of box k	*/
				ff_temp[k][i]		= ff[k][i];
				vv_temp[k][i]		= vv[k][i];
			}
							
				/*======================================*/
				/* Tranfer data from box k to temp	*/
				/*======================================*/
				en_temp[k]			= en[k];


			MPI_Recv(atom[k],   box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank+1,3,MPI_COMM_WORLD,&status);
			MPI_Recv(atnopbc[k],box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank+1,4,MPI_COMM_WORLD,&status);
			MPI_Recv(ff[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank+1,5,MPI_COMM_WORLD,&status);
			MPI_Recv(vv[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank+1,6,MPI_COMM_WORLD,&status);
			MPI_Recv(&en[k],     sizeof(struct energy),MPI_CHAR,mpi.my_rank+1,11,MPI_COMM_WORLD,&status);

			uu[k]			= vv[k];
			ordparam[k].d_nc= swap_accept.d_nc;
			
			int new_bin_1	= (int) ((ordparam[k].d_nc - sim_dos[k].l_begin)/sim_dos[k].l_width);

			MPI_Send(atom_temp[k],   box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank+1, 7,MPI_COMM_WORLD);
			MPI_Send(atnopbc_temp[k],box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank+1, 8,MPI_COMM_WORLD);
			MPI_Send(ff_temp[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank+1, 9,MPI_COMM_WORLD);
			MPI_Send(vv_temp[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank+1,10,MPI_COMM_WORLD);
			MPI_Send(&en_temp[k],     sizeof(struct energy),MPI_CHAR,mpi.my_rank+1,12,MPI_COMM_WORLD);
			
		

			dos_hist[k][new_bin_1].h_of_l ++;
			dos_hist[k][new_bin_1].g_of_l += sim_dos[k].mod_f;
			flat_histogram(k);
			#ifdef NLIST
				 nblist(k);									
			#endif
			if(sim_dos[k].mod_f <MERGE_F && ((iter%(MERGE_N*(unsigned long)sim_hyb.cyc_swap)==0) || 
						(iter%(((unsigned long)(MERGE_N+1))*(unsigned long)sim_hyb.cyc_swap)==0))){
				 MPI_Send(dos_hist[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank+1,14,MPI_COMM_WORLD);
				 MPI_Recv(dos_hist[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank+1,15,MPI_COMM_WORLD,&status);
			}
			end_to_end(k,SITE1,SITE2);
			return 1;
		}// if swap accepted then for the sending box
	 }// end of mpi.rank !=p-1
	 
	 else if(mpi.my_rank!=0 && flag==0)
	 {
		/* Recv data from my_rank-1 */
		struct msg_swap_request swap_request;
		struct msg_swap_accept swap_accept;
		/* Recv data from my_rank-1 */
		MPI_Recv(&swap_request,sizeof(swap_request),MPI_CHAR,mpi.my_rank-1,1,MPI_COMM_WORLD,&status);
		double l2 = swap_request.d_nc;
		double mod_f = swap_request.mod_f;

		/* --------------------------------------------	*/
		/*	This onwards added later for on the fly	*/
		/*	merging of g_of_e			*/
		/* -------------------------------------------- */
		double l_begin	=	swap_request.l_begin;
		double l_end	=	swap_request.l_end;

		if ((l2 >= sim_dos[k].l_begin && l2 < sim_dos[k].l_end) && (l1 >= swap_request.l_begin && l1 < swap_request.l_end)) 
		{
			swap_accept.accept = 1;
			swap_accept.d_nc = l1;
			MPI_Send(&swap_accept,sizeof(swap_accept),   MPI_CHAR,mpi.my_rank-1,2,MPI_COMM_WORLD);
			MPI_Send(atom[k],   box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank-1,3,MPI_COMM_WORLD);
			MPI_Send(atnopbc[k],box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank-1,4,MPI_COMM_WORLD);
			MPI_Send(ff[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank-1,5,MPI_COMM_WORLD);
			MPI_Send(vv[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank-1,6,MPI_COMM_WORLD);
			MPI_Send(&en[k],     sizeof(struct energy),MPI_CHAR,mpi.my_rank-1,11,MPI_COMM_WORLD);
			ordparam[k].d_nc	= l2;
			int new_bin_1	= (int) ((ordparam[k].d_nc - sim_dos[k].l_begin)/sim_dos[k].l_width);
			MPI_Recv(atom[k],   box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank-1, 7,MPI_COMM_WORLD,&status);
			MPI_Recv(atnopbc[k],box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank-1, 8,MPI_COMM_WORLD,&status);
			MPI_Recv(ff[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank-1, 9,MPI_COMM_WORLD,&status);
			MPI_Recv(vv[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank-1,10,MPI_COMM_WORLD,&status);
			MPI_Recv(&en[k],     sizeof(struct energy),MPI_CHAR,mpi.my_rank-1,12,MPI_COMM_WORLD,&status);
			uu[k]	=	vv[k];
			dos_hist[k][new_bin_1].h_of_l ++;
			dos_hist[k][new_bin_1].g_of_l += sim_dos[k].mod_f;
			flat_histogram(k);
			#ifdef NLIST
					 nblist(k);									
			#endif
			
			if(mod_f<MERGE_F && ((iter%(MERGE_N*(unsigned long)sim_hyb.cyc_swap)==0) || 
						(iter%(((unsigned long)(MERGE_N+1))*(unsigned long)sim_hyb.cyc_swap)==0))){
				MPI_Recv(dos_hist_temp[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank-1,14,MPI_COMM_WORLD,&status);
				if (sim_dos[k].mod_f <MERGE_F) merge_gofl(k,l_begin,l_end, mod_f);
				MPI_Send(dos_hist_temp[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank-1,15,MPI_COMM_WORLD);
			}
			end_to_end(k,SITE1,SITE2);
			return 1;
		}//end of accept for receiving box
		else {
			swap_accept.accept=0;
			swap_accept.d_nc=l2;// l2 belongs to myrank-1
			MPI_Send(&swap_accept,sizeof(swap_accept),   MPI_CHAR,mpi.my_rank-1,2,MPI_COMM_WORLD);

			dos_hist[k][old_bin_1].h_of_l ++;
			dos_hist[k][old_bin_1].g_of_l += sim_dos[k].mod_f;
			flat_histogram(k);
			
			if(mod_f<MERGE_F && ((iter%(MERGE_N*(unsigned long)sim_hyb.cyc_swap)==0) || 
						(iter%(((unsigned long)(MERGE_N+1))*(unsigned long)sim_hyb.cyc_swap)==0))){
				MPI_Recv(dos_hist_temp[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank-1,14,MPI_COMM_WORLD,&status);
				if (sim_dos[k].mod_f <MERGE_F) merge_gofl(k,l_begin,l_end,mod_f);
				MPI_Send(dos_hist_temp[k],NBINS*sizeof(struct hist),MPI_CHAR,mpi.my_rank-1,15,MPI_COMM_WORLD);
			}
			end_to_end(k,SITE1,SITE2);
			return 0;
		}// end of reject for receiving box
	 }// end of rank!=0 loop
	 return 2;
}

void merge_gofl(int k, double l_begin, double l_end, double mod_f){
	int l_bins		= (int)((l_end - l_begin)/sim_dos[k].l_width+0.5);// Number of bins in my_rank -1 histogram
	int o_bins		= (int)((l_end - sim_dos[k].l_begin)/sim_dos[k].l_width);

	double sum=0.0;
	for(int i=0; i<o_bins; i++){
		double diff = dos_hist_temp[k][l_bins - o_bins +i].g_of_l - dos_hist[k][i].g_of_l;
		sum			+= diff;
	}
	sum /= o_bins;
	for(int j=0; j< sim_dos[k].l_bins; j++){
		dos_hist[k][j].g_of_l += sum;
	}
	for (int i=0; i<o_bins; i++){
		dos_hist[k][i].g_of_l = (dos_hist[k][i].g_of_l + dos_hist_temp[k][l_bins - o_bins +i].g_of_l)*0.5;
		dos_hist_temp[k][l_bins - o_bins +i].g_of_l = dos_hist[k][i].g_of_l;
	}

	
}
#endif//XEDOS
#endif//MPI
