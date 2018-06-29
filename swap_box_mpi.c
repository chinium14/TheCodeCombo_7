#ifdef MPI
/* ========================================================================= */
/* swap_box_mpi.cpp							  */
/*									   */
/* This subroutine is called from repexch.cpp to swap the conformations      */
/* in adjacent boxes. It is the parallel compliment of swap_box.	     */
/* Differnt processes are executed on different processors.		  */
/* ========================================================================= */

#include "defines.h"

void nblist	 (int);
double ran2	 (void);
#ifdef STATUS
void curr_status    (int,int);
#endif
int swap_box_mpi(int flag)
{
	int k = 0;
	/*------------------------------------------------------*/
	/* The energy of this process is saved in e1 whereas	*/
	/* The energy of the partner process is saved in e2	*/
	/*------------------------------------------------------*/
	MPI_Status status;
	#ifdef STATUS
	curr_status(k,5);
	#endif
	if(mpi.my_rank!=mpi.p-1 && flag ==1)//flag=1 mean the box is sending a request to rank+1
	{
		struct msg_swap_request swap_request;
            #ifndef SASAREX
		swap_request.potens  =       en[k].potens;
            #else
		swap_request.potens  =      en[k].esasab - en[k].esasa;  // esasab: the same configuration in a new solvent condition
                                                                         // esasa: the corrent 
            #endif
		swap_request.kT	     =       sim.kT[k];
		struct msg_swap_accept swap_accept;

		/* ------------------------------------------------ */
		/* If I am not the last process and if my flag == 1 */
		/* then send data/swap request to my_rank +1	    */
		/* ------------------------------------------------ */
		MPI_Send(&swap_request,sizeof(swap_request),MPI_CHAR,mpi.my_rank+1,1,MPI_COMM_WORLD);

		/* ------------------------------------------------ */
		/* Receive the information whether move is accepted */
		/* or rejected.						*/
		/* ------------------------------------------------ */
		MPI_Recv(&swap_accept, sizeof(swap_accept), MPI_CHAR,mpi.my_rank+1,2,MPI_COMM_WORLD,&status);
		/* ------------------------------------------------ */
		/* If move is rejected no more communication needed */
		/* ------------------------------------------------ */
		if(swap_accept.accept==0)
		{
			return 0;
		}
		/* ------------------------------------------------ */
		/* If move is accepted then create a back up of	    */
		/* of current structures and then receive the data	*/
		/* from the partner process and after that send the */
		/* back up info to the partner to complete the swap */
		/* ------------------------------------------------ */
		else if (swap_accept.accept==1)
		{

			/*======================================*/
			/* Tranfer data from box k to temp			*/
			/*======================================*/
			for(int i =0; i< box[k].boxns; i++)
			{
				atom_temp[k][i]		= atom[k][i];	      /* Back up coordinates of box k			*/
				atnopbc_temp[k][i]	= atnopbc[k][i];  /* Back up coordinates of box k			*/
				ff_temp[k][i]		= ff[k][i];
				vv_temp[k][i]		= vv[k][i];
			}
			en_temp[k]			= en[k];
			#ifdef PRESSURE
			pvir_temp[k]		      = pvir[k];
			#endif
			/*======================================*/
			/* Receive data from other box	  */
			/*======================================*/
			MPI_Recv(atom[k],   box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank+1,3,MPI_COMM_WORLD,&status);
			MPI_Recv(atnopbc[k],box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank+1,4,MPI_COMM_WORLD,&status);
			MPI_Recv(ff[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank+1,5,MPI_COMM_WORLD,&status);
			MPI_Recv(vv[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank+1,6,MPI_COMM_WORLD,&status);
			MPI_Recv(&en[k],    sizeof(struct energy),MPI_CHAR,mpi.my_rank+1,11,MPI_COMM_WORLD,&status);
			#ifdef PRESSURE
			MPI_Recv(&pvir[k],    sizeof(struct virial),MPI_CHAR,mpi.my_rank+1,13,MPI_COMM_WORLD,&status);
			#endif
			for(int i=0; i< box[k].boxns; i++)
			{
				vv[k][i].x  /= swap_accept.scale;
				vv[k][i].y  /= swap_accept.scale;
				vv[k][i].z  /= swap_accept.scale;
				uu[k][i]     = vv[k][i];
			}
			// en[k].potens	= swap_accept.potens;
			/*======================================*/
			/* Send data to other box	       */
			/*======================================*/
			MPI_Send(atom_temp[k],   box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank+1, 7,MPI_COMM_WORLD);
			MPI_Send(atnopbc_temp[k],box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank+1, 8,MPI_COMM_WORLD);
			MPI_Send(ff_temp[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank+1, 9,MPI_COMM_WORLD);
			MPI_Send(vv_temp[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank+1,10,MPI_COMM_WORLD);
			MPI_Send(&en_temp[k],     sizeof(struct energy),MPI_CHAR,mpi.my_rank+1,12,MPI_COMM_WORLD);
			#ifdef PRESSURE
			//MPI_Send(&pvir[k],    sizeof(struct virial),MPI_CHAR,mpi.my_rank+1,14,MPI_COMM_WORLD,&status);
			MPI_Send(&pvir[k],    sizeof(struct virial),MPI_CHAR,mpi.my_rank+1,14,MPI_COMM_WORLD);
			#endif
			#ifdef NLIST
			nblist(k);
			#endif
			return 1;
		}// if swap accepted then for the sending box
	}// end of mpi.rank !=p-1
	else if(mpi.my_rank!=0 && flag==0)//flag 0 means the box is receiving a request from rank-1
	{
		struct msg_swap_request swap_request;
		struct msg_swap_accept swap_accept;
		/*==========================*/
		/* Recv data from my_rank-1 */
		/*==========================*/
		MPI_Recv(&swap_request,sizeof(swap_request),MPI_CHAR,mpi.my_rank-1,1,MPI_COMM_WORLD,&status);
		/* ---------------------------------------- */
		/* Accept criterion for the swap	    */
		/* ---------------------------------------- */
		double e2 = swap_request.potens;
                #ifndef SASAREX
		   double delta_energy = e2 - en[k].potens;
   		   double T2 = swap_request.kT;
		   double delta_beta = 1.0/T2-1.0/ sim.kT[k];
		   double accep_crit = exp(delta_beta*delta_energy);
                #else
                   double delta_energy = -(e2 + en[k].esasab - en[k].esasa);
                   double accep_crit = exp((1.0/sim.kT[k])*delta_energy);
                #endif
		/* ---------------------------------------- */
		/* Evaluate if swap will take place.	*/
		/* ---------------------------------------- */
		if(accep_crit > ran2())
		{//Swap is accepted
			swap_accept.accept = 1;
                #ifndef SASAREX
			swap_accept.scale  = sqrt(sim.kT[k]/T2);
                #else
                        swap_accept.scale  = 1;
                #endif
			/* ====================== */
			/* Send this box info to  */
			/* other box.	     */
			/* ====================== */
			MPI_Send(&swap_accept,sizeof(swap_accept),   MPI_CHAR,mpi.my_rank-1,2,MPI_COMM_WORLD);
			MPI_Send(atom[k],   box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank-1,3,MPI_COMM_WORLD);
			MPI_Send(atnopbc[k],box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank-1,4,MPI_COMM_WORLD);
			MPI_Send(ff[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank-1,5,MPI_COMM_WORLD);
			MPI_Send(vv[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank-1,6,MPI_COMM_WORLD);
			MPI_Send(&en[k],     sizeof(struct energy),MPI_CHAR,mpi.my_rank-1,11,MPI_COMM_WORLD);
			#ifdef PRESSURE
			//MPI_Send(&pvir[k],    sizeof(struct virial),MPI_CHAR,mpi.my_rank-1,13,MPI_COMM_WORLD,&status);
			MPI_Send(&pvir[k],    sizeof(struct virial),MPI_CHAR,mpi.my_rank-1,13,MPI_COMM_WORLD);
			#endif
			/* ====================== */
			/* Recieve new box info   */
			/* from other box for     */
			/* this box.	      */
			/* ====================== */
			MPI_Recv(atom[k],   box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank-1, 7,MPI_COMM_WORLD,&status);
			MPI_Recv(atnopbc[k],box[k].boxns * sizeof(struct atoms),MPI_CHAR,mpi.my_rank-1, 8,MPI_COMM_WORLD,&status);
			MPI_Recv(ff[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank-1, 9,MPI_COMM_WORLD,&status);
			MPI_Recv(vv[k],     box[k].boxns * sizeof(struct veloc),MPI_CHAR,mpi.my_rank-1,10,MPI_COMM_WORLD,&status);
			MPI_Recv(&en[k],     sizeof(struct energy),MPI_CHAR,mpi.my_rank-1,12,MPI_COMM_WORLD,&status);
			#ifdef PRESSURE
			MPI_Recv(&pvir[k],    sizeof(struct virial),MPI_CHAR,mpi.my_rank-1,14,MPI_COMM_WORLD,&status);
			#endif

			/* ===================== */
			/* Scale the velocities  */
			/* as prescribe by the   */
			/* repexch algorithm.  	*/
			/* ===================== */
			for(int i=0; i< box[k].boxns; i++)
			{
				vv[k][i].x  *= swap_accept.scale;
				vv[k][i].y  *= swap_accept.scale;
				vv[k][i].z  *= swap_accept.scale;
				uu[k][i]     = vv[k][i];
			}
			#ifdef NLIST
			nblist(k);						
			#endif
			return 1;
		}//end of accept for receiving box
		else
		{ //swap is rejected
			swap_accept.accept=0;
			swap_accept.potens=e2;// e2 belongs to myrank-1
			MPI_Send(&swap_accept,sizeof(swap_accept),   MPI_CHAR,mpi.my_rank-1,2,MPI_COMM_WORLD);
			return 0;
		}// end of reject for receiving box
	}// end of rank!=0 loop
	return 2;
}

#endif
