/* ======================================================================== */
/* clist.cpp                                                                */
/* Written by Thomas Knotts 26 Apr 2003                                     */
/*                                                                          */
/*		This file contains subroutines that initialize a cell list and      */
/* determine the neighbors of cells.  The grid is set up so cells have      */
/* a length of about CELL_LENGTH anstroms.  This method currently only      */
/* works with square simulation boxes that do not change shape (ie.. you    */
/* cannot use PR_NPT with this option.  It also will assign the same number */
/* of cells                                                                 */
/* ======================================================================== */
#ifdef CLIST
#include "defines.h"


int icell(int,int,int,int,int,int);

/* ************************************************************************ */
/* ======================== Begin Subroutine ============================== */
/* ************************************************************************ */
void init_cell(void){

	/* ================================================================== */
	/* Allocate memory for the cell neighbor list.                        */
	/* ================================================================== */
	clistn	= (struct clistn_struct*) calloc(sim.NB,sizeof(struct clistn_struct));
		if(clistn == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for clistn\n"); exit(11);}
	head	= (int**) calloc(sim.NB,sizeof(int*));
		if(head == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for head\n"); exit(11);}
	clist	= (int**) calloc(sim.NB,sizeof(int*)); 
		if(clist == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for clist\n"); exit(11);}
	cell_no	= (int**) calloc(sim.NB,sizeof(int*));
		if(cell_no == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for cell_no\n"); exit(11);}
	for(int k=0; k<sim.NB; k++){
		cell_no[k]	= (int*) calloc(box[k].boxns,sizeof(int));
			if(cell_no[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for cell_no[k]\n"); exit(11);}
	}
	for(int k=0; k<sim.NB; k++){	
	/* ================================================================== */
	/* Determine the number of cell per side(M) so that the length of     */
	/* the cell is about CELL_LENGTH anstrom.  Then, determine the        */
	/* total number of boxes, the number of shells of neighboring cells,  */
	/* and allocate the necessary memory.                                 */
	/* ================================================================== */
		int Mx, My, Mz;

		Mx = (int)(box[k].lx/CELL_LENGTH_X);
		My = (int)(box[k].ly/CELL_LENGTH_Y);
		Mz = (int)(box[k].lz/CELL_LENGTH_Z);
		clistn[k].Mx = Mx;
		clistn[k].My = My;
		clistn[k].Mz = Mz;
		int ncell = Mx * My * Mz;
		clistn[k].ncell = ncell;
		double cell_lengthx = box[k].lx/(double)(Mx);
		double cell_lengthy = box[k].ly/(double)(My);
		double cell_lengthz = box[k].lz/(double)(Mz);

		double least_cell_length = cell_lengthx;
		if(least_cell_length > cell_lengthy) least_cell_length=cell_lengthy;
		if(least_cell_length > cell_lengthz) least_cell_length=cell_lengthz;

	/* ------------------------------------------------------------------ */
	/* Determine the number of shells of boxes to go out to check for     */
	/* neighbors.                                                         */
	/* ------------------------------------------------------------------ */
		int nshellx, nshelly, nshellz;

#ifdef EWALD
		double rc = sim.rclist_ew;
#else
		double rc = sim.rclist;
#endif

		if(cell_lengthx < rc){
			nshellx = (int)(rc/cell_lengthx);
			double test = rc/cell_lengthx-(double)nshellx;
			if(test > 0.001) nshellx++;
		}
		else
			nshellx = 1;

		if(cell_lengthy < rc){
			nshelly = (int)(rc/cell_lengthy);
			double test = rc/cell_lengthy-(double)nshelly;
			if(test > 0.001) nshelly++;
		}
		else
			nshelly = 1;

		if(cell_lengthz < rc){
			nshellz = (int)(rc/cell_lengthz);
			double test = rc/cell_lengthz-(double)nshellz;
			if(test > 0.001) nshellz++;
		}
		else
			nshellz = 1;

	/* ================================================================== */
	/* Allocate memory                                                    */
	/* ================================================================== */
		
		head[k]			= (int*) calloc(ncell+MEXTRA,sizeof(int*));
			if(head[k]  == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for head[k]\n"); exit(11);}
		clist[k]		= (int*) calloc(ncell+box[k].boxns+MEXTRA,sizeof(int*)); 
			if(clist[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for clist[k]\n"); exit(11);}
		clistn[k].count	= (int*) calloc(ncell+MEXTRA,sizeof(int));
			if(clistn[k].count  == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for clistn[k].count\n"); exit(11);}
		clistn[k].list	= (int**) calloc(ncell+MEXTRA,sizeof(int*));
			if(clistn[k].list  == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for clistn[k].list\n"); exit(11);}
	/* ------------------------------------------------ */
	/* Check to make sure the shell length is at least  */
	/* as large as the cut off for the neighborlist,    */
	/* but not too large to overcount.                  */
	/* ------------------------------------------------ */
		double shell_lengthx = (double)(nshellx)*cell_lengthx;
		if(cell_lengthx < rc){
			if(shell_lengthx<rc){
				printf("Problem, shell length x=%f is less than rcut=%f!!\n", shell_lengthx, rc);
				exit(840);
			}
		}
		double shell_lengthy = (double)(nshelly)*cell_lengthy;
		if(cell_lengthy < rc){
			if(shell_lengthy<rc){
				printf("Problem, shell length y=%f is less than rcut=%f!!\n", shell_lengthy, rc);
				exit(840);
			}
		}
		double shell_lengthz = (double)(nshellz)*cell_lengthz;
		if(cell_lengthz < rc){
			if(shell_lengthz<rc){
				printf("Problem, shell length z=%f is less than rcut=%f!!\n", shell_lengthz, rc);
				exit(840);
			}
		}

		if(nshellx+1 > Mx/2){
			printf("Problem, you must make CELL_LENGTH_X (defines.h) smaller or you will overcount atoms!!\n");
			exit(840);
		}
		if(nshelly+1 > My/2){
			printf("Problem, you must make CELL_LENGTH_Y (defines.h) smaller or you will overcount atoms!!\n");
			exit(840);
		}
		if(nshellz+1 > Mz/2){
			printf("Problem, you must make CELL_LENGTH_Z (defines.h) smaller or you will overcount atoms!!\n");
			exit(840);
		}

	/* ================================================================== */
	/* Allocate the memory needed for the cell neighbor list.             */
	/* ================================================================== */
		int size = (int)((2*nshellx+1)*(2*nshelly+1)*(2*nshellz+1));
		if(size ==-1) size = 1;
		for(int i=0; i<ncell; i++){
			clistn[k].list[i] = (int*)calloc(size,sizeof(int));
				if(clistn[k].list[i]  == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for clistn[k].list[i]\n"); exit(11);}
		}
	/* ================================================================== */
	/* Loop to determine the neighbors of each box.                       */
	/* ================================================================== */
		for(int iz=0; iz<Mz; iz++){
			for(int iy=0; iy<My; iy++){
				for(int ix=0; ix<Mx; ix++){
					int cell = icell(ix,iy,iz,Mx,My,Mz);
					for(int dz=-nshellz; dz<=nshellz; dz++){		
						for(int dy=-nshelly; dy<=nshelly; dy++){
							for(int dx=-nshellx; dx<=nshellx; dx++){
								//if(abs(dx)+abs(dy)+abs(dz)!=0){
									int cell_nabor = icell(ix+dx,iy+dy,iz+dz,Mx,My,Mz);
									//if(cell_nabor > cell){
										clistn[k].list[cell][clistn[k].count[cell]]=cell_nabor;
										clistn[k].count[cell]++;
									//}//if(cell_nabor > cell)

								//}//if(dx*dy*dz!=0)
							}
						}
					}
				}
			}
		}
	}//for k
}//end subroutine init_cell

int icell(int ix, int iy, int iz, int Mx, int My, int Mz){
	if(ix<0) ix += Mx;
	if(iy<0) iy += My;
	if(iz<0) iz += Mz;
	return(ix%Mx + iy%My * Mx + iz%Mz * Mx * My);
}

#endif
