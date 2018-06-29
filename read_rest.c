/* ======================================================================== */
/* read_rest.cpp	                                                        */
/* This suboutine reads the information about the atoms to be restrained    */
/* if REST is defined.  It allocates the memory and creates the list of     */
/* restrained atoms to be used in crestraint.cpp.                           */ 
/*																		    */
/* Written by Thomas Knotts 27 Aug 03                                       */
/* ======================================================================== */
#ifdef REST
#include "defines.h"

void read_rest(void){
  
  char tt[150],file_name[50];
  FILE *file_ptr;
  int sa;

  int n_total;



  /* ====================================================== */
  /* Allocate the memory for the rest1 and rest2            */
  /* global variables.                                      */
  /* ====================================================== */
  rest1 = (struct rest_struct_1*)calloc(sim.NB,sizeof(struct rest_struct_1));
    if(rest1 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest1\n"); exit(11);}
  rest2 = (struct rest_struct_2*)calloc(sim.NB,sizeof(struct rest_struct_2));
    if(rest2 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest2\n"); exit(11);}

  /* ====================================================== */
  /* Open the rest#.input file, read the number of          */
  /* restrained sites.  Allocate a temporary array for      */
  /* one atom restaints from two atom restraints.  Then     */
  /* read the file into this structure.                     */
  /* ====================================================== */


  for(int k=0; k<sim.NB; k++){

  struct{
	int n1;			// number of one atom restraints
	int n2;			// number of two atom restraints
	int *site1;		// site # 1 of restraint
	int *site2;		// site # 2 of restraint
	double *k;		
	double *req;	// equilibirum distance for two atom restraints
	double *pos_x;
	double *pos_y;
	double *pos_z;
  } rest_temp;

	  
	  sprintf(file_name,"./INPUT/rest%d.input",k);
	  if( NULL == (file_ptr=fopen(file_name,"r")) ) {
		  fprintf(stdout,"ERROR: Input file %s does not exist!!!\n",file_name);
		  exit(1);
	  }
  
	  fgets(tt,150,file_ptr);
      fgets(tt,150,file_ptr);

	  fscanf(file_ptr,"%d", &n_total);			fgets(tt,150,file_ptr);


	  rest_temp.n1 = 0;
	  rest_temp.n2 = 0;
	  rest_temp.site1	= (int*)calloc(n_total+1,sizeof(int));
	    if(rest_temp.site1 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest_temp.site1\n"); exit(11);}
	  rest_temp.site2 = (int*)calloc(n_total+1,sizeof(int));
	    if(rest_temp.site2 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest_temp.site2\n"); exit(11);}
	  rest_temp.k		= (double*)calloc(n_total+1,sizeof(double));
	    if(rest_temp.k == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest_temp.k\n"); exit(11);}
	  rest_temp.req		= (double*)calloc(n_total+1,sizeof(double));
	    if(rest_temp.req == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest_temp.req\n"); exit(11);}
	  rest_temp.pos_x		= (double*)calloc(n_total+1,sizeof(double));
/* The following added so that coords for restrain can be read from rest.inp file	*/

	    if(rest_temp.pos_x == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest_temp.pos_x\n"); exit(11);}
	  rest_temp.pos_y		= (double*)calloc(n_total+1,sizeof(double));
	    if(rest_temp.pos_y == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest_temp.pos_y\n"); exit(11);}
	  rest_temp.pos_z		= (double*)calloc(n_total+1,sizeof(double));
	    if(rest_temp.pos_z == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest_temp.pos_z\n"); exit(11);}

	  for(int i=0; i<n_total; i++){
		  fscanf(file_ptr,"%d",&sa);
		  fscanf(file_ptr,"%d",&rest_temp.site1[i]);
		  fscanf(file_ptr,"%d",&rest_temp.site2[i]);
		  fscanf(file_ptr,"%lf",&rest_temp.k[i]);
		  fscanf(file_ptr,"%lf",&rest_temp.req[i]);
		  if(rest_temp.site2[i]==0){
			  fscanf(file_ptr,"%lf",&rest_temp.pos_x[i]);
			  fscanf(file_ptr,"%lf",&rest_temp.pos_y[i]);
			  fscanf(file_ptr,"%lf",&rest_temp.pos_z[i]);			  
		  }
		  fgets(tt,150,file_ptr);
		  if(rest_temp.site1[i] == rest_temp.site2[i]) rest_temp.n1++;
		  else if(rest_temp.site2[i]==0) rest_temp.n1++;
		  else rest_temp.n2++;

		  rest_temp.site1[i]--;
		  rest_temp.site2[i]--;
	  }

  /* ====================================================== */
  /* Allocate the memory for the rest1 and rest2            */
  /* global variables.                                      */
  /* ====================================================== */
	  rest1[k].site		=(int*)calloc(rest_temp.n1+1,sizeof(int));
	    if(rest1[k].site == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest1[k].site\n"); exit(11);}
	  rest1[k].k		=(double*)calloc(rest_temp.n1+1,sizeof(double));
	    if(rest1[k].k == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest1[k].k\n"); exit(11);}
	  rest1[k].x0		=(double*)calloc(rest_temp.n1+1,sizeof(double));
	    if(rest1[k].x0 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest1[k].x0\n"); exit(11);}
	  rest1[k].y0		=(double*)calloc(rest_temp.n1+1,sizeof(double));
	    if(rest1[k].y0 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest1[k].y0\n"); exit(11);}
	  rest1[k].z0		=(double*)calloc(rest_temp.n1+1,sizeof(double));
	    if(rest1[k].z0 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest1[k].z0\n"); exit(11);}
      rest1[k].req		=(double*)calloc(rest_temp.n1+1,sizeof(double));
	    if(rest1[k].req == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest1[k].req\n"); exit(11);}


	  rest2[k].site1	=(int*)calloc(rest_temp.n2+1,sizeof(int));
	    if(rest2[k].site1 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest2[k].site1\n"); exit(11);}
	  rest2[k].site2	=(int*)calloc(rest_temp.n2+1,sizeof(int));
	    if(rest2[k].site2 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest2[k].site2\n"); exit(11);}
	  rest2[k].k		=(double*)calloc(rest_temp.n2+1,sizeof(double));
	    if(rest2[k].k == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest2[k].k\n"); exit(11);}
	  rest2[k].req		=(double*)calloc(rest_temp.n2+1,sizeof(double));
	    if(rest2[k].req == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rest2[k].req\n"); exit(11);}

  /* ====================================================== */
  /* Now, loop through the rest_temp structure, check to    */
  /* see if the restraint is with one or two atoms, and     */
  /* fill the rest1 and rest2 global arrays.                */      
  /* ====================================================== */
	  int n1=0;
	  int n2=0;

	  rest1[k].n = rest_temp.n1;
	  rest2[k].n = rest_temp.n2;

	  for(int i=0; i<n_total; i++){
		  if(rest_temp.site2[i] == -1){
			  rest1[k].site[n1] = rest_temp.site1[i];
			  rest1[k].k[n1] = rest_temp.k[i] * 4.184;
			  rest1[k].x0[n1] = rest_temp.pos_x[i];
			  rest1[k].y0[n1] = rest_temp.pos_y[i];
			  rest1[k].z0[n1] = rest_temp.pos_z[i];
			  rest1[k].req[n1]	 = rest_temp.req[i];
			  n1++;
		  }
		  else if(rest_temp.site1[i] == rest_temp.site2[i]){
			  rest1[k].site[n1] = rest_temp.site1[i];
			  rest1[k].k[n1] = rest_temp.k[i] * 4.184;
			  rest1[k].x0[n1] = atom[k][rest1[k].site[n1]].x;
			  rest1[k].y0[n1] = atom[k][rest1[k].site[n1]].y;
			  rest1[k].z0[n1] = atom[k][rest1[k].site[n1]].z;
			  rest1[k].req[n1]	 = rest_temp.req[i];
			  n1++;
		  }
		  else{
			  rest2[k].site1[n2] = rest_temp.site1[i];
			  rest2[k].site2[n2] = rest_temp.site2[i];
			  rest2[k].k[n2]	 = rest_temp.k[i] * 4.184;
			  rest2[k].req[n2]	 = rest_temp.req[i];
			  n2++;
		  }
	  }
	  if((n1+n2)!= n_total){
		  printf("There is an error with the restrained sites: n_total (%i) does not equal n1+n2 (%i+%i=%i).", n_total, n1, n2, n1+n2);
		  exit(501);
	  }
	  if(n1 != rest_temp.n1){
		  printf("There is an error with the restrained sites(n1).");
		  exit(502);
	  }
	  if(n2 != rest_temp.n2){
		  printf("There is an error with the restrained sites(n2).");
		  exit(503);
	  }
	  fclose(file_ptr);
  }//k
}
#endif


