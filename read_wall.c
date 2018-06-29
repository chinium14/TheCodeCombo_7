#ifdef WALL
/* ======================================================================== */
/* read_wall.cpp															*/
/* This suboutine reads the information for walls in the system and			*/
/* does the initialization necessary for these interactions.                */
/*																			*/
/* Written by Thomas Knotts 11 Mar 04                                       */
/* ======================================================================== */
#include "defines.h"

void read_wall(void){
  
  wall    = (struct wall_struct*) calloc(sim.NB,sizeof(struct wall_struct));
		if(wall == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for wall\n"); exit(11);}


  char tt[150],file_name[50];
  FILE *file_ptr;
  int boxnum;  
  /* ====================================================== */
  /* Open the wall.input file and read in parameters.       */
  /* ====================================================== */
  sprintf(file_name,"./INPUT/wall.input");
  if( NULL == (file_ptr=fopen(file_name,"r")) ) {
    fprintf(stdout,"ERROR: Input file %s does not exist!!!\n",file_name);
    exit(1);
  }
    
  fgets(tt,150,file_ptr);
  fgets(tt,150,file_ptr);

  for(int k=0; k<sim.NB; k++){    

	  fscanf(file_ptr,"%d",&boxnum);		fgets(tt,150,file_ptr);
	  if(boxnum != k){ fprintf(stdout, "ERROR: Box numbers (%i) in golik.input are not correct! (should be %i)\n", boxnum, k); exit(010);}

	  fscanf(file_ptr,"%d",&wall[k].n);		fgets(tt,150,file_ptr);
	  if(wall[k].n > 1){
		  printf("Multiple walls is not supported right now. (In XEDOS the check to see if a particle is past the wall is only valid for one wall.\n");
		  exit(34);
	  }
	  wall[k].z = (double*)calloc(wall[k].n,sizeof(double));
		 if(wall[k].z == NULL){fprintf(stdout, "ERROR: cannot allocate memory for wall[k].z\n"); exit(11);}
	  wall[k].nsites = (int*)calloc(wall[k].n,sizeof(int));
		 if(wall[k].nsites == NULL){fprintf(stdout, "ERROR: cannot allocate memory for wall[k].nsites\n"); exit(11);}
	  wall[k].num_density = (double**)calloc(wall[k].n,sizeof(double*));
		 if(wall[k].num_density == NULL){fprintf(stdout, "ERROR: cannot allocate memory for wall[k].num_density\n"); exit(11);}
	  wall[k].eps = (double**)calloc(wall[k].n,sizeof(double*));
		 if(wall[k].eps == NULL){fprintf(stdout, "ERROR: cannot allocate memory for wall[k].eps\n"); exit(11);}
	  wall[k].sig = (double**)calloc(wall[k].n,sizeof(double*));
		 if(wall[k].sig == NULL){fprintf(stdout, "ERROR: cannot allocate memory for wall[k].sig\n"); exit(11);}
	  wall[k].hydro_index = (double**)calloc(wall[k].n,sizeof(double*));
		 if(wall[k].hydro_index == NULL){fprintf(stdout, "ERROR: cannot allocate memory for wall[k].hydro_index\n"); exit(11);}

	  for(int i=0; i<wall[k].n; i++){
		fscanf(file_ptr,"%lf",&wall[k].z[i]);                fgets(tt,150,file_ptr);
		fscanf(file_ptr,"%d",&wall[k].nsites[i]);           fgets(tt,150,file_ptr);

		wall[k].num_density[i] = (double*)calloc(wall[k].nsites[i],sizeof(double));
		  if(wall[k].num_density[i] == NULL){fprintf(stdout, "ERROR: cannot allocate memory for wall[k].num_density[i]\n"); exit(11);}
		wall[k].eps[i] = (double*)calloc(wall[k].nsites[i],sizeof(double));
		  if(wall[k].eps[i] == NULL){fprintf(stdout, "ERROR: cannot allocate memory for wall[k].eps[i]\n"); exit(11);}
		wall[k].sig[i] = (double*)calloc(wall[k].nsites[i],sizeof(double));
		  if(wall[k].sig[i] == NULL){fprintf(stdout, "ERROR: cannot allocate memory for wall[k].sig[i]\n"); exit(11);}
		wall[k].hydro_index[i] = (double*)calloc(wall[k].nsites[i],sizeof(double));
		  if(wall[k].hydro_index[i] == NULL){fprintf(stdout, "ERROR: cannot allocate memory for wall[k].hydro_index[i]\n"); exit(11);}


		for(int j=0; j<wall[k].nsites[i]; j++){
		  fscanf(file_ptr,"%lf",&wall[k].num_density[i][j]);
		  fscanf(file_ptr,"%lf",&wall[k].sig[i][j]);
		  fscanf(file_ptr,"%lf",&wall[k].eps[i][j]);
		  wall[k].eps[i][j] *= 4.184;
		  fscanf(file_ptr,"%lf",&wall[k].hydro_index[i][j]);
		  fgets(tt,150,file_ptr);	
		}
	  }
  
		fscanf(file_ptr, "%d %lf %lf",&wall[k].angle_site,&wall[k].angle,&wall[k].angle_k); fgets(tt,150,file_ptr);
		wall[k].angle  *= PI/180.0;
		wall[k].angle_k *= 4.184;
//		printf("PI is %lf \n", PI);
//		printf("wall[k].angle is %lf\n", wall[k].angle);
//		printf("wall[k].angle_k is %lf\n", wall[k].angle_k);
		wall[k].angle_site--;
//		printf("wall[k].angle_site is %d\n", wall[k].angle_site);
		//printf("angle site = %d\n",wall[k].angle_site);
/* This was an attempt to get the wall mixing rules into ljset.  It might work, but cwall wasn't coded to use this.
   In the future, if you want to use this, you need to change golik[k].d_repulsive below to the site-specific sigma.
   Then, you need to change cwall to use ljset[k].pott[][].eps/sig.  You will have to address STYPE/no STYPE.
		for(int a=0; a<wall[k].n; a++){
			for(int n=0;  n < wall[k].nsites[a]; n++){
				for(int b=0; b<box[k].boxns; b++){
					ljset[k].pot[box[k].boxns+a][b].eps = sqrt(pott[k][b].eps * wall[k].eps[a][n]);
					ljset[k].pot[box[k].boxns+a][b].sig = 0.5*(golik[k].d_repulsive + wall[k].sig[a][n]);
				}
			}
		}*/

  }

	  fclose(file_ptr);

}
#endif
