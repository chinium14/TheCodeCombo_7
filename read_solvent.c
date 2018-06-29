/* ======================================================================== */
/* read_rest.cpp	                                                        */
/* This suboutine reads the information about the atoms to be restrained    */
/* if REST is defined.  It allocates the memory and creates the list of     */
/* restrained atoms to be used in crestraint.cpp.                           */ 
/*																		    */
/* Written by Thomas Knotts 27 Aug 03                                       */
/* ======================================================================== */
#ifdef SASA
#ifdef BGO
#include "defines.h"

void read_solvent(void){
  
  char tt[150],file_name[50];
  FILE *file_ptr;
#ifdef SASAREX
  double concentration_local;
#endif 
  for(int k=0; k<sim.NB; k++){
#ifdef MPI
          sprintf(file_name,"./INPUT/solvent%d.input",mpi.my_rank);
#else	  
	  sprintf(file_name,"./INPUT/solvent%d.input",k);
#endif
	  if( NULL == (file_ptr=fopen(file_name,"r")) ) {
		  fprintf(stdout,"ERROR: Input file %s does not exist!!!\n",file_name);
		  exit(1);
	  }
  
	  fgets(tt,150,file_ptr);
          fgets(tt,150,file_ptr);

          fscanf(file_ptr,"%d", &solvent_type);                          fgets(tt,150,file_ptr);
//	  fscanf(file_ptr,"%lf", &concentration[k]);			fgets(tt,150,file_ptr);
#ifdef SASAREX
	  fscanf(file_ptr,"%lf", &concentration_local);			fgets(tt,150,file_ptr);
#else
	  fscanf(file_ptr,"%lf", &concentration[k]);			fgets(tt,150,file_ptr);
          // printf("machine # %d: concentration [%d] is %lf\n", mpi.my_rank, k, concentration[k]); 
#endif
	  fclose(file_ptr);
  }//k
#ifdef SASAREX
  MPI_Status status;
  int root = 0;
  if(mpi.my_rank == root){ 
          for (int k=1; k<mpi.p; k++){
              MPI_Recv(&concentration[k],1,MPI_DOUBLE,k,1,MPI_COMM_WORLD,&status);
          }     
  }
  else{
          MPI_Send(&concentration_local,1,MPI_DOUBLE,root,1,MPI_COMM_WORLD);
  }    
  MPI_Bcast(concentration,24,MPI_DOUBLE,root,MPI_COMM_WORLD);    
 // printf("machine # %d: concentration [%d] is %lf\n", mpi.my_rank, mpi.my_rank, concentration[mpi.my_rank]); 
#endif
#ifdef MEMINSIDE
  char in_list[50];
  FILE *in_ptr;
          sprintf(in_list,"./INPUT/in_list.inp");
          if( NULL == (in_ptr=fopen(in_list,"r")) ) {
                  fprintf(stdout,"ERROR: Input file %s does not exist!!!\n",in_list);
                  exit(1);
          }
          fscanf(in_ptr,"%d",&number_in); fgets(tt,150,in_ptr);
          mem_in_list = (int*) calloc(number_in, sizeof(int));
          for (int i=0; i<number_in; i++){
             fscanf(in_ptr,"%d", &mem_in_list[i]); fgets(tt,150,in_ptr);
          }
          fclose(in_ptr);
#endif
}
#endif
#endif
