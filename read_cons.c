/* ======================================================================== */
/* read_cons.cpp	                                                        */
/* This suboutine reads the information about the atoms to be constrained   */
/* if KONS is defined.  It allocates the memory and creates the list of    */
/* constrained atoms to be used in forces.cpp, force_long.cpp, and          */
/* force_short.cpp to set the velocities and forces to zero.                */ 
/*																		    */
/* Written by Thomas Knotts 26 Aug 03                                       */
/* ======================================================================== */
#ifdef KONS
#include "defines.h"

void read_cons(void){
  
  char tt[150],file_name[50];
  FILE *file_ptr;
  int sa;

  /* ====================================================== */
  /* Open the cons#.input file, read the number of          */
  /* constrained sites, allocated memory, and generate the  */
  /* cons[k].site[i] array.                                 */
  /* ====================================================== */

  cons = (struct cons_struct*) calloc(sim.NB,sizeof(struct cons_struct));
  if(cons == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for cons\n"); exit(11);}

  for(int k=0; k<sim.NB; k++){
	  sprintf(file_name,"./INPUT/cons%d.input",k);
	  if( NULL == (file_ptr=fopen(file_name,"r")) ) {
		  fprintf(stdout,"ERROR: Input file %s does not exist!!!\n",file_name);
		  exit(1);
	  }
  

	  fgets(tt,150,file_ptr);
      fgets(tt,150,file_ptr);

      fscanf(file_ptr,"%d", &cons[k].n);			fgets(tt,150,file_ptr);

      cons[k].site = (int*) calloc(cons[k].n,sizeof(int));
	  if(cons[k].site == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for cons[k].site\n"); exit(11);}

      for(int i=0; i<cons[k].n; i++){
		  fscanf(file_ptr,"%d", &sa);
	      fscanf(file_ptr,"%d", &cons[k].site[i]);			fgets(tt,150,file_ptr);
		  cons[k].site[i] -= 1;
	  }
	  fclose(file_ptr);
  }
}
#endif

