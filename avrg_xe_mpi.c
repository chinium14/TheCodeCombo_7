/* ========================================================================= */
/* avrg_xe_mpi.cpp                                                           */
/*                                                                           */
/* ========================================================================= */

#include "defines.h"
#ifdef MPI
#ifdef XEDOS
void nblist		  (int);
double ran2		  (void);
void flat_histogram(int);
void average_gofl (int);
void end_to_end	  (int,int,int);


int avrg_xe_mpi(int iter, int flag, int ibox)
{
	int k =ibox;
	MPI_Status status;
	int send_rank = mpi.my_rank;
	int recv_rank = mpi.p - mpi.my_rank-1;

	if(flag ==0){
		MPI_Send(dos_hist[k],NBINS*sizeof(struct hist),MPI_CHAR,recv_rank,17,MPI_COMM_WORLD);
    	MPI_Recv(dos_hist[k],NBINS*sizeof(struct hist),MPI_CHAR,recv_rank,18,MPI_COMM_WORLD,&status);
		return 1;
	 }
	else if(flag==1){
		MPI_Recv(dos_hist_temp[k],NBINS*sizeof(struct hist),MPI_CHAR,recv_rank,17,MPI_COMM_WORLD,&status);
		average_gofl(k);
		MPI_Send(dos_hist_temp[k],NBINS*sizeof(struct hist),MPI_CHAR,recv_rank,18,MPI_COMM_WORLD);
		return 1;
	}
	else return 0;
}

void average_gofl(int k){
	int l_bins = sim_dos[k].l_bins;
	for (int i =0; i<l_bins; i++){
		dos_hist[k][i].g_of_l = (dos_hist[k][i].g_of_l + dos_hist_temp[k][l_bins -i-1].g_of_l)*0.5;
		dos_hist[k][i].h_of_l = (int)((dos_hist[k][i].h_of_l + dos_hist_temp[k][l_bins -i-1].h_of_l)*0.5);
	}
	for (int i =0; i<l_bins; i++){
		dos_hist_temp[k][i].g_of_l = dos_hist[k][l_bins-i-1].g_of_l;
		dos_hist_temp[k][i].h_of_l = dos_hist[k][l_bins-i-1].h_of_l;
	}

}
#endif//XEDOS
#endif//MPI
