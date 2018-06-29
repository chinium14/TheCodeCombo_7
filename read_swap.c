#ifdef TDXEDOS
/* ======================================================================== */
/* read_swap.cpp		                                                    */
/* This suboutine reads the swap#.cpp file and calculates the               */
/* probabilities of swapping.  It is used when doing 2DXEDOS.               */
/*                                                                          */
/* Written by Thomas Knotts 4 June 04                                       */	
/* ======================================================================== */
#include "defines.h"
void read_swap(void){
  int temp1;
  char tt[150], name[100];
  FILE *file_ptr;

  /* ======================================================== */
  /* Determine the number of boxes/processors and allocate    */
  /* memory.                                                  */
  /* ======================================================== */
  #ifdef MPI
    swapping.n_boxes = mpi.p;
  #else
    swapping.n_boxes = sim.NB;
  #endif
if(swapping.n_boxes > 1){
  swapping.n_nabors    = (int*)calloc(swapping.n_boxes,sizeof(int));
    if(swapping.n_nabors == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for swapping.n_nabors\n"); exit(11);}
  swapping.nabor       = (int**)calloc(swapping.n_boxes,sizeof(int*));
    if(swapping.nabor == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for swapping.nabor\n"); exit(11);}
  for(int i=0; i<swapping.n_boxes; i++){
    swapping.nabor[i]  = (int*)calloc(8,sizeof(int));
	  if(swapping.nabor[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for swapping.nabor[i]\n"); exit(11);}
  }
  swapping.pick_box_prob = (double*)calloc(swapping.n_boxes,sizeof(double));
    if(swapping.pick_box_prob == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for swapping.pick_box_prob\n"); exit(11);}
  /* ======================================================== */
  /* Open the file and read the information.                  */
  /* ======================================================== */

  sprintf(name, "./INPUT/swap.input");


  file_ptr = fopen(name,"r");
  if( NULL == file_ptr){
    printf("Input file 2Dxedos.input does not exist!!!\n");
    exit(1);
  }
  fgets(tt,150,file_ptr);
  fgets(tt,150,file_ptr);
  fgets(tt,150,file_ptr);
  fgets(tt,150,file_ptr);


  for(int k=0; k<swapping.n_boxes; k++){
    fscanf(file_ptr,"%d",&temp1);                      fgets(tt,150,file_ptr);
	if(temp1 != k){ printf("The box number (%i) is not correct in swap.input (should be %i)", temp1, k); exit(10);}
	fscanf(file_ptr,"%d",&swapping.n_nabors[k]);       fgets(tt,150,file_ptr);
    for(int i=0; i<swapping.n_nabors[k]; i++){
	  fscanf(file_ptr,"%d",&swapping.nabor[k][i]);    fgets(tt,150,file_ptr);
	}
	fgets(tt,150,file_ptr);
  }

  /* ======================================================== */
  /* Determine the probability with which to pick each box    */
  /* in proposing a swap and create the probability table.    */
  /* ======================================================== */
  int total_nabors = 0;
  for(int k=0; k<swapping.n_boxes; k++) total_nabors += swapping.n_nabors[k];

  swapping.pick_box_prob[0] = (double)swapping.n_nabors[0]/(double)total_nabors;
  for(int k=1; k<swapping.n_boxes; k++){
    swapping.pick_box_prob[k] = swapping.pick_box_prob[k-1]+(double)swapping.n_nabors[k]/(double)total_nabors;
  }
  fclose(file_ptr);
}
}
#endif

