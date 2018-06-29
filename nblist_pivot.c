/* ======================================================================== */
/* nblist_pivot.cpp                                                               */
/*                                                                          */
/*		This subroutine sets up the neighbor-list for each molecule. It is  */
/* called in main only if NLIST is defined.                                 */
/*                                                                          */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

#ifdef NLIST
#ifdef CLIST
int icell2(double,double,double,double,double,double,int,int,int);
#endif


/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void nblist_pivot (int ibox, int ub)
{
  int k = ibox;

  /* ================================================================== */
  /*                                                                    */
  /* Get and assign information for minimum image convention and        */
  /* neighborlist cutoff.                                               */
  /*                                                                    */
  /* ================================================================== */
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;

  double rcut = sim.rclist;
  double rcut2 = sim.rclist * sim.rclist;
  int ns = box[k].boxns;

#ifdef SLIST
  double rcslist	=	sim.rc + 0.5*(rcut - sim.rc);
  double rcslist2	=	rcslist*rcslist;
#endif
#ifdef EWALD
  double rcut_ew = sim.rclist_ew;
  double rcut2_ew = rcut_ew * rcut_ew;
#endif  
#ifdef SASA
  double rcutsasa2 = 9.0*9.0;
#endif

  /* ================================================================== */
  /* Zero out count accumulators.                                       */
  /* ================================================================== */
   for(int i=0; i<ub; i++) {
	  nlist[k].count[i] = 0;

#ifdef SASA
	  nlist[k].count_sasa[i] =0;
#endif
#ifdef EWALD
	  nlist_ew[k].count[i] = 0;
#endif
  }
#ifdef SLIST
   int nss = molec.lsite[molec.Nsolute-1];
      for(int i=0; i<nss; i++) {
		slist[k].count[i] = 0;
	  }
#endif
  /* ================================================================== */
  /* Loops to create neighbor list.                                     */
  /* ================================================================== */
#ifndef CLIST
   
   for(int a=0; a<ub; a++) {

  /* ----------------------------------------------- */
  /* Assign position of particle a.                  */
  /* ----------------------------------------------- */    
    double ax = atom[k][a].x;
    double ay = atom[k][a].y;
    double az = atom[k][a].z;

  /* ----------------------------------------------- */
  /* Loop around atoms "above" particle a in array.  */
  /* ----------------------------------------------- */
    for(int b=a+1; b<ns; b++) {

  /* ----------------------------------------------- */
  /* Calculate the x,y,z distances between           */
  /* particles a & b and apply minimum image.        */
  /* ----------------------------------------------- */
      double dx = fabs(ax - atom[k][b].x);
		if(dx >  hx) dx -= lx;
		if(dx > rcut) continue;	  
	  double dy = fabs(ay - atom[k][b].y);
		if(dy >  hy) dy -= ly;
		if(dy > rcut) continue;
      double dz = fabs(az - atom[k][b].z);
		if(dz >  hz) dz -= lz;
		if(dz > rcut) continue;

  /* ----------------------------------------------- */
  /* Calculate the distance between particles a & b. */
  /* If the distance is less than the cutoff,        */
  /* add the particle to the list. If it is greater, */
  /* go to the next atom in the loop.                */
  /* ----------------------------------------------- */
      double dr2 = dx*dx + dy*dy + dz*dz;

#ifdef EWALD
	  if(dr2 < rcut2_ew){
		nlist_ew[k].list[a][nlist_ew[k].count[a]] = b;
#ifdef STYPE
		nlist_ew[k].xflag[a][nlist_ew[k].count[a]] = 0;  //the flag is zero if the pair is not part of a bond, bend, or torsion
  /* ----------------------------------------------- */
  /* Search through the list of excluded             */
  /* interactions and assign the appropriate flag    */
  /* to the neighboring pair.                        */
  /* ----------------------------------------------- */
	if(xintN[k]){
		int index = (a*b) % NPRIME;
		int n = table[k][index].n;
		for(int r=0; r<n; r++){
			int ra = table[k][index].a[r];
			int rb = table[k][index].b[r];
			if(((a == ra) && (b == rb)) || ((a == rb) && (b == ra))){
				nlist_ew[k].xflag[a][nlist_ew[k].count[a]] = xint[k][table[k][index].p[r]].xflag; 
			}//if(((a == ra) && (b == rb)) || ((a == rb) && (b == ra)))
		
		}//for r


	}//
#endif
	
	
	if(dr2 < rcut2){
		nlist[k].list[a][nlist[k].count[a]] = b;
#ifdef STYPE
		nlist[k].xflag[a][nlist[k].count[a]] = nlist_ew[k].xflag[a][nlist_ew[k].count[a]];
#endif
		nlist[k].count[a]++;
		if(nlist[k].count[a] > max_num_neighbors){
			printf("The neighborlist array is too small.  Increase max_num_neighbors to at least %i in allocate.cpp or check your density!\n", nlist[k].count[a]);
			exit(5);
		}
	}

	nlist_ew[k].count[a]++;
	if(nlist_ew[k].count[a] > max_num_neighbors){
		printf("The ewald neighborlist array is too small.  Increase max_num_neighbors to at least %i in allocate.cpp or check your density!\n", nlist_ew[k].count[a]);
		exit(5);
	}


#else //ifdef EWALD

      if(dr2 < rcut2) {

      
	nlist[k].list[a][nlist[k].count[a]] = b;

#ifdef STYPE
	nlist[k].xflag[a][nlist[k].count[a]] = 0;  //the flag is zero if the pair is not part of a bond, bend, or torsion
  /* ----------------------------------------------- */
  /* Search through the list of excluded             */
  /* interactions and assign the appropriate flag    */
  /* to the neighboring pair.                        */
  /* ----------------------------------------------- */
	if(xintN[k]){
		int index = (a*b) % NPRIME;
		int n = table[k][index].n;
		for(int r=0; r<n; r++){
			int ra = table[k][index].a[r];
			int rb = table[k][index].b[r];
			if(((a == ra) && (b == rb)) || ((a == rb) && (b == ra))){
				nlist[k].xflag[a][nlist[k].count[a]] = xint[k][table[k][index].p[r]].xflag; 
#ifdef DNA_GOLIK
        nlist[k].sig[a][nlist[k].count[a]]=xint[k][table[k][index].p[r]].sig;
#endif
#ifdef SASA
				if (dr2 < rcutsasa2){
					nlist[k].xflag_sasa[a][nlist[k].count_sasa[a]] = xint[k][table[k][index].p[r]].xflag;
				}
#endif
			}//if(((a == ra) && (b == rb)) || ((a == rb) && (b == ra)))
		
		}//for r


	}//
#endif
	nlist[k].count[a]++;
	if(nlist[k].count[a] > max_num_neighbors){
		printf("The neighborlist array is too small.  Increase max_num_neighbors to at least %i in allocate.cpp or check your density!\n", nlist[k].count[a]);
		exit(5);
	}

#endif //ewald
  /* ----------------------------------------------- */
  /* If SASA is defined, there are two neighbor      */
  /* lists.  One for nonbonded interactions, and one */
  /* to calculate the solvent accessible surface     */
  /* area.                                           */
  /* ----------------------------------------------- */
#ifdef SASA
	if (dr2 < rcutsasa2){
		nlist[k].list_sasa[a][nlist[k].count_sasa[a]] = b;
		nlist[k].count_sasa[a]++;
	}
#endif
#ifdef SLIST
	if (dr2 < rcslist2){
		if((a<molec.lsite[molec.Nsolute-1]) && (strcmp(atom[k][b].name,"OH2")) == 0){
			slist[k].list[a][slist[k].count[a]]=b;
			slist[k].count[a] ++;
		}
	}
#endif

	
      }// if(dr2 < rcut2)

    }//end b loop
  }// end a loop


#endif //ifndef CLIST

#ifdef CLIST
  /* ================================================================== */
  /* Loop to assign sites to cells.                                     */
  /* ================================================================== */
  int ncell = clistn[k].ncell;
  int Mx = clistn[k].Mx;
  int My = clistn[k].My;
  int Mz = clistn[k].Mz;
  int cell;

  for(int icell=0; icell<ncell; icell++) head[k][icell]=-1;
  for(int ilist=0; ilist<ncell+ns; ilist++) clist[k][ilist] = -1;
	
  for(int i=0; i<ns; i++){
	cell = icell2((atom[k][i].x + hx), (atom[k][i].y + hy),(atom[k][i].z + hz),lx,ly,lz,Mx,My,Mz);
	cell_no[k][i] = cell;
	clist[k][i] = head[k][cell];
	head[k][cell] = i;
  }
		  
  /* ================================================================== */
  /* Assign the neighbor list from the linked list.                     */
  /* ================================================================== */
  /* ------------------------------------------------------------ */
  /* Loop over N-1 particles.                                     */
  /* ------------------------------------------------------------ */
	for(int a=0; a<ub; a++) {
		double ax = atom[k][a].x;
		double ay = atom[k][a].y;
		double az = atom[k][a].z;
  /* ------------------------------------------------------------ */
  /* Loop over the neighbors in the current cell which are below  */
  /* particle a.                                                  */
  /* ------------------------------------------------------------ */
		cell = cell_no[k][a];
		
		for(int cnabor=0; cnabor<clistn[k].count[cell]; cnabor++){
			int b = head[k][clistn[k].list[cell][cnabor]];
			while(b!=-1){
				if(b>a){
					double dx = fabs(ax - atom[k][b].x);
						if(dx >  hx) dx -= lx;
						if(dx > rcut){
							b = clist[k][b];
							continue;	  
						}
					double dy = fabs(ay - atom[k][b].y);
						if(dy >  hy) dy -= ly;
						if(dy > rcut){
							b = clist[k][b];
							continue;	  
						}
					double dz = fabs(az - atom[k][b].z);
						if(dz >  hz) dz -= lz;
						if(dz > rcut){
							b = clist[k][b];
							continue;	  
						}


  /* ----------------------------------------------- */
  /* Calculate the distance between particles a & b. */
  /* If the distance is less than the cutoff,        */
  /* add the particle to the list. If it is greater, */
  /* go to the next atom in the loop.                */
  /* ----------------------------------------------- */
					double dr2 = dx*dx + dy*dy + dz*dz;
#ifdef EWALD
	  if(dr2 < rcut2_ew){
		nlist_ew[k].list[a][nlist_ew[k].count[a]] = b;
#ifdef STYPE
		nlist_ew[k].xflag[a][nlist_ew[k].count[a]] = 0;  //the flag is zero if the pair is not part of a bond, bend, or torsion
  /* ----------------------------------------------- */
  /* Search through the list of excluded             */
  /* interactions and assign the appropriate flag    */
  /* to the neighboring pair.                        */
  /* ----------------------------------------------- */
	if(xintN[k]){
		int index = (a*b) % NPRIME;
		int n = table[k][index].n;
		for(int r=0; r<n; r++){
			int ra = table[k][index].a[r];
			int rb = table[k][index].b[r];
			if(((a == ra) && (b == rb)) || ((a == rb) && (b == ra))){
				nlist_ew[k].xflag[a][nlist_ew[k].count[a]] = xint[k][table[k][index].p[r]].xflag; 
			}//if(((a == ra) && (b == rb)) || ((a == rb) && (b == ra)))
		
		}//for r


	}//
#endif
	
	
	if(dr2 < rcut2){
		nlist[k].list[a][nlist[k].count[a]] = b;
#ifdef STYPE
		nlist[k].xflag[a][nlist[k].count[a]] = nlist_ew[k].xflag[a][nlist_ew[k].count[a]];
#endif
		nlist[k].count[a]++;
		if(nlist[k].count[a] > max_num_neighbors){
			printf("The neighborlist array is too small.  Increase max_num_neighbors to at least %i in allocate.cpp or check your density!\n", nlist[k].count[a]);
			exit(5);
		}

	}

	nlist_ew[k].count[a]++;
	if(nlist_ew[k].count[a] > max_num_neighbors){
		printf("The ewald neighborlist array is too small.  Increase max_num_neighbors to at least %i in allocate.cpp or check your density!\n", nlist_ew[k].count[a]);
		exit(5);
	}

#else //ifdef EWALD

					if(dr2 < rcut2) {
						nlist[k].list[a][nlist[k].count[a]] = b;
#ifdef STYPE
						nlist[k].xflag[a][nlist[k].count[a]] = 0;  //the flag is zero if the pair is not part of a bond, bend, or torsion
  /* ----------------------------------------------- */
  /* Search through the list of excluded             */
  /* interactions and assign the appropriate flag    */
  /* to the neighboring pair.                        */
  /* ----------------------------------------------- */
						if(xintN[k]){
							int index = (a*b) % NPRIME;
							int n = table[k][index].n;
							for(int r=0; r<n; r++){
								int ra = table[k][index].a[r];
								int rb = table[k][index].b[r];
								if(((a == ra) && (b == rb)) || ((a == rb) && (b == ra))){
									nlist[k].xflag[a][nlist[k].count[a]] = xint[k][table[k][index].p[r]].xflag; 
#ifdef DNA_GOLIK
                  nlist[k].sig[a][nlist[k].count[a]]=xint[k][table[k][index].p[r]].sig;
#endif
#ifdef SASA
									if (dr2 < rcutsasa2){
										nlist[k].xflag_sasa[a][nlist[k].count_sasa[a]] = xint[k][table[k][index].p[r]].xflag;
									}
#endif
								}//if(((a == ra) && (b == rb)) || ((a == rb) && (b == ra)))
							}//for r
						}//xintN
#endif
						nlist[k].count[a]++;
						if(nlist[k].count[a] > max_num_neighbors){
							printf("The neighborlist array is too small.  Increase max_num_neighbors to at least %i in allocate.cpp or check your density!\n", nlist[k].count[a]);
							exit(5);
						}

#endif //EWALD
  /* ----------------------------------------------- */
  /* If SASA is defined, there are two neighbor      */
  /* lists.  One for nonbonded interactions, and one */
  /* to calculate the solvent accessible surface     */
  /* area.                                           */
  /* ----------------------------------------------- */
#ifdef SASA
						if (dr2 < rcutsasa2){
							nlist[k].list_sasa[a][nlist[k].count_sasa[a]] = b;
							nlist[k].count_sasa[a]++;
						}
#endif
#ifdef SLIST
	if (dr2 < rcslist2){
		if((a<molec.lsite[molec.Nsolute-1]) && (strcmp(atom[k][b].name,"OH2")) == 0){
			slist[k].list[a][slist[k].count[a]]=b;
			slist[k].count[a] ++;
		}
	}
#endif
					}// if(dr2 < rcut2)
					b = clist[k][b];
				}//if(b<0a)
				else b = clist[k][b];
			}//while b!=-1
		}//for cnabor
	}//for a

#endif //CLIST
#ifdef DLIST
	for(int i=0; i<ns; i++){
		rx0[k][i] = atom[k][i].x;
		ry0[k][i] = atom[k][i].y;
		rz0[k][i] = atom[k][i].z;
	}
	nl_flag[k] = 0;
#endif
}

#endif //NLIST 
