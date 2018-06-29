#if defined WWALL || defined SPHERE
#include "defines.h"
/* ======================================================================== */
/* ======================================================================== */
void set_go_params(int ibox);
int is_bond(int,int,int);
void hydro_index_bgo(void){
//  FILE * hi;
//  hi=fopen("hydro_index.dat","w");
  double hydro_index[47]={0.0,-0.4,1.8,-0.8,2.5,4.2,-0.7,4.5,1.6,1.9,-3.5,-3.5,3.8,-3.9,-3.5,
    -3.5,-4.5,-3.2,2.8,-1.3,-0.9,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  bgo    = (struct bgo_struct*) calloc(sim.NB,sizeof(struct bgo_struct));
                if(bgo == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for bgo\n"); exit(11);}
  int boxnum;
  for(int k=0; k<sim.NB; k++){
          bgo[k].hydro_index  = (double*)calloc(box[k].boxns,sizeof(double));
                if(bgo[k].hydro_index == NULL){fprintf(stdout, "ERROR: cannot allocate memory for bgo[k].hydro_index\n"); exit(11);}

          int rescount = 0;
          int atcount = 0;
          for(int m=0; m<sim.NC; m++) {
                for(int i=0; i<bp[k][m].nbox; i++) {
                  for (int kk=0; kk<mol[m].Nres; kk++) {
                        for(int j=0; j<residue[k][rescount].Nsite; j++) {
		                  bgo[k].hydro_index[atcount] = hydro_index[residue[k][rescount].type];
//				  fprintf(hi,"hydro_index is %f\n",bgo[k].hydro_index[atcount]);
                  		  atcount++;
              		}//for j
                	rescount++;
                  }//for kk
                }//for i
          }//for m
  }//for k
//  fclose(hi);
}
#endif
