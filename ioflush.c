/* ======================================================================== */
/* ioflush.cpp                                                              */
/*                                                                          */
/*       This subroutine checks to see if the buffers containing the output */
/* data are full and if they are writes them to file.                       */
/*                Parameters:	flag: If it is 0, check before flushing     */
/*									  buffer, if it is anything else, flush */
/*                                    without checking.                     */
/* Written by Thomas Knotts, 8 Sep 04                                       */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/* ======================== Begin Subroutine ============================== */
/* ************************************************************************ */

void ioflush(int flag) 
{
	FILE *file_ptr;
	char file_name[100];
	int k_rank;

	for(int k=0; k<sim.NB; k++){

		#ifdef MPI
		k_rank = mpi.my_rank;
		#else
		k_rank = k;
		#endif

		/* =========================================== */
		/* Write out if longer than 1 hour since last  */
		/* output.                                     */
		/* =========================================== */
		time_t now;
		time(&now);
		double max_io_time = 10;	//max time between io in seconds


		
		if(flag == 0){

			if(simul_ptr[k] - simul_buff[k] > OUTBUFF || difftime(now,simul_time[k])> max_io_time){
				sprintf(file_name,"./OUTPUT/BOX%d/simul%d.output",k_rank,k_rank);
				file_ptr = fopen(file_name,"a");
				fwrite(simul_buff[k],sizeof(char),(simul_ptr[k] - simul_buff[k]),file_ptr);
				fclose(file_ptr);
				simul_ptr[k] = simul_buff[k];
				time(&simul_time[k]);
			}

			if(ener_ptr[k] - ener_buff[k] > OUTBUFF || difftime(now,ener_time[k])> max_io_time){
				sprintf(file_name,"./OUTPUT/BOX%d/ener_box%d.output",k_rank,k_rank);
				file_ptr = fopen(file_name,"a");
				fwrite(ener_buff[k],sizeof(char),(ener_ptr[k] - ener_buff[k]),file_ptr);
				fclose(file_ptr);
				ener_ptr[k] = ener_buff[k];
				time(&ener_time[k]);
			}

			if(ord_ptr[k] - ord_buff[k] > OUTBUFF || difftime(now,ord_time[k])> max_io_time){
				sprintf(file_name,"./OUTPUT/BOX%d/order%d.output",k_rank,k_rank);
				file_ptr = fopen(file_name,"a");
				fwrite(ord_buff[k],sizeof(char),(ord_ptr[k] - ord_buff[k]),file_ptr);
				fclose(file_ptr);
				ord_ptr[k] = ord_buff[k];
				time(&ord_time[k]);
			}

			if(swap_ptr[k] - swap_buff[k] > OUTBUFF || difftime(now,swap_time[k])> max_io_time){

/*				#if defined(DOS)   || defined(MMDOS)  || defined(CTDOS)   || \
					defined(XEDOS) || defined(FXEDOS) || defined(TDXEDOS)
					#ifdef MPI
					sprintf(file_name,"./OUTPUT/DOS/swap_dos%d.out",k_rank);
					#else
					sprintf(file_name,"./OUTPUT/DOS/swap_dos.out");
					#endif
				#else
				sprintf(file_name,"./OUTPUT/swap_box.out");
				#endif*/
        sprintf(file_name,"./OUTPUT/swap%d.out",k_rank);

				if((sim.ID == 2) || (sim.ID==4) || (sim.ID == 5) || (sim.ID == 7) || (sim.ID == 9) || (sim.ID == 10) 
					|| (sim.ID == 11) || (sim.ID == 12)){
					file_ptr = fopen(file_name,"a");
					fwrite(swap_buff[k],sizeof(char),(swap_ptr[k] - swap_buff[k]),file_ptr);
					fclose(file_ptr);
					swap_ptr[k] = swap_buff[k];
					time(&swap_time[k]);
				}
			}

			#if defined(SMD) || defined(NSMD)
			if(smd_ptr[k] - smd_buff[k] > OUTBUFF || difftime(now,smd_time[k])> max_io_time){
				sprintf(file_name,"./OUTPUT/BOX%d/smd%d.output",k_rank,k_rank);
				file_ptr = fopen(file_name,"a");
				fwrite(smd_buff[k],sizeof(char),(smd_ptr[k] - smd_buff[k]),file_ptr);
				fclose(file_ptr);
				smd_ptr[k] = smd_buff[k];
				time(&smd_time[k]);
			}
			#endif
		
	  
			
		}
		else{

			sprintf(file_name,"./OUTPUT/BOX%d/simul%d.output",k_rank,k_rank);
			file_ptr = fopen(file_name,"a");
			fwrite(simul_buff[k],sizeof(char),(simul_ptr[k] - simul_buff[k]),file_ptr);
			fclose(file_ptr);
			simul_ptr[k] = simul_buff[k];
			time(&simul_time[k]);

			sprintf(file_name,"./OUTPUT/BOX%d/ener_box%d.output",k_rank,k_rank);
			file_ptr = fopen(file_name,"a");
			fwrite(ener_buff[k],sizeof(char),(ener_ptr[k] - ener_buff[k]),file_ptr);
			fclose(file_ptr);
			ener_ptr[k] = ener_buff[k];
			time(&ener_time[k]);

			sprintf(file_name,"./OUTPUT/BOX%d/order%d.output",k_rank,k_rank);
			file_ptr = fopen(file_name,"a");
			fwrite(ord_buff[k],sizeof(char),(ord_ptr[k] - ord_buff[k]),file_ptr);
			fclose(file_ptr);
			ord_ptr[k] = ord_buff[k];
			time(&ord_time[k]);

			#if defined(SMD) || defined(NSMD)
			sprintf(file_name,"./OUTPUT/BOX%d/smd%d.output",k_rank,k_rank);
			file_ptr = fopen(file_name,"a");
			fwrite(smd_buff[k],sizeof(char),(smd_ptr[k] - smd_buff[k]),file_ptr);
			fclose(file_ptr);
			smd_ptr[k] = smd_buff[k];
			time(&smd_time[k]);
			#endif


			/*#if defined(DOS)   || defined(MMDOS)  || defined(CTDOS)   || \
				defined(XEDOS) || defined(FXEDOS) || defined(TDXEDOS)
                #ifdef MPI
               sprintf(file_name,"./OUTPUT/swap%d.out",k_rank);
                #else
                sprintf(file_name,"./OUTPUT/DOS/swap_dos.out");
                #endif			
			#else
			sprintf(file_name,"./OUTPUT/swap_box.out");
			#endif*/

			if((sim.ID == 2) || (sim.ID == 4) || (sim.ID == 5) || (sim.ID == 7) || (sim.ID == 9) || (sim.ID == 10) 
				|| (sim.ID == 11) || (sim.ID == 12)){
				file_ptr = fopen(file_name,"a");
				fwrite(swap_buff[k],sizeof(char),(swap_ptr[k] - swap_buff[k]),file_ptr);
				fclose(file_ptr);
				swap_ptr[k] = swap_buff[k];
				time(&swap_time[k]);
			}

			
		}//if(flag == 0)
	}//k
	
}

