/* ======================================================================== */
/* setlj.cpp                                                                */
/*                                                                          */
/*		This subroutine sets the Lennard-Jones parameters for each          */
/* interaction. It applies the Lorentz-Bertheloth mixing rules and stores   */
/* the parameters for the interaction of particles types a and b in         */
/* ljset[a][b].  This the the structure used in cnonbond and                */
/* cononbond_nblist to calculate the LJ interactions.  The hash table,      */
/* used for the excluded interactions, is also set up here.                 */
/*                                                                          */
/* Passed Parameters:                                                       */
/*                                                                          */
/* ======================================================================== */
#include "defines.h"


/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void setlj (void)
{
#ifndef STYPE
  for(int k=0; k<sim.NB; k++){
	  /* ================================================================== */
	  /* Loop around particle interactions and apply the mixing rules.      */
	  /* ================================================================== */
    for(int a=0; a<box[k].boxns; a++){
      #ifndef CLIST
      for(int b=a; b<box[k].boxns; b++){
      #else
      for(int b=0; b<box[k].boxns; b++){
      #endif

        #ifdef EWALD
        ljset[k].pot[a][b].bond_flag = 0;
        ljset[k].pot[a][b].bend_flag = 0;
        ljset[k].pot[a][b].tors_flag = 0;
        ljset[k].pot[a][b].exNB_flag = 0;
        #endif		
        double eps = sqrt(pott[k][a].eps * pott[k][b].eps);
        double sig = 0.5 * (pott[k][a].sig + pott[k][b].sig);
	double qq  = pott[k][a].qq * pott[k][b].qq;
	//int id = 0; //0 means there are no speacial nbfix, native contacts, etc.
	#ifdef SASA
        #ifndef BGO
	ljset[k].pot[a][b].pab = 0.3516;// parameter for non bonded pairs for implicit solvent model
        #endif
        #ifdef BGO
	ljset[k].pot[a][b].pab = 0.5872;// parameter for non bonded pairs for BGO solvent model
        #endif
	#endif
        /* ================================================================== */
        /* If a and b are bonded, remove the LJ interaction.                  */
        /* ================================================================== */
        for (int n=0; n<bondN[k]; n++){
          if ((a == bond[k][n].a && b== bond[k][n].b) || (a == bond[k][n].b && b== bond[k][n].a)){
            eps = 0.0;
            sig = 0.0;
            qq  = 0.0;
	    #ifdef EWALD
            ljset[k].pot[a][b].bond_flag = 1;
	    #endif
            #ifdef SASA
            #ifndef BGO
            ljset[k].pot[a][b].pab = 0.8875;
            #endif
            #ifdef BGO
            ljset[k].pot[a][b].pab = 1.24310;   
            #endif
            #endif
            break;
          }
        }
        /* ================================================================== */
        /* If a and b are part of a bend, remove the LJ interaction.          */
        /* ================================================================== */
        for (int n=0; n<bendN[k]; n++) {
          if ((a == bend[k][n].a && b== bend[k][n].c) || (a == bend[k][n].c && b== bend[k][n].a)){
            eps = 0.0;
            sig = 0.0;
            qq  = 0.0;
	    #ifdef EWALD
            ljset[k].pot[a][b].bend_flag = 1;
	    #endif
            #ifdef SASA
            #ifdef BGO
            ljset[k].pot[a][b].pab = -0.15956;
            #endif
            #endif
            break;
          }
        }

        #ifdef S14
        /* ================================================================== */
        /* If a and b are part of a 1-4 interaction, use the special LJ       */
        /* parameters and scale qq by E14FAC.                                 */
        /* ================================================================== */
        for (int n=0; n<in14N[k]; n++) {
          if ((a == in14[k][n].a && b== in14[k][n].d) || (a == in14[k][n].d && b== in14[k][n].a)){
            eps = sqrt(pott14[k][a].eps * pott14[k][b].eps); 
            sig = 0.5 * (pott14[k][a].sig + pott14[k][b].sig);
            qq *= E14FAC;
	    #ifdef EWALD
            ljset[k].pot[a][b].tors_flag = 1;
	    #endif
            #ifdef SASA
            #ifdef BGO
            ljset[k].pot[a][b].pab = 0.65562;
            #endif
            #endif
            break;
          }
        }
        #endif //S14

        /* ================================================================== */
        /* If a and b are part of a NBFIX, use the specified interaction      */
        /* parameters.                                                        */
        /* ================================================================== */
        for (int n=0; n<NUM_NBFIX; n++){
	  int ta = atnopbc[k][a].atomid;
	  int tb = atnopbc[k][b].atomid;
          if ((ta == nbfix_prop[n].type1 && tb == nbfix_prop[n].type2) ||
	      (tb == nbfix_prop[n].type1 && ta == nbfix_prop[n].type2)){
            eps = nbfix_prop[n].eps*-4.184;
            sig = nbfix_prop[n].sig;
            //id = 10;
            break;
          }
        }
        /* ================================================================== */
        /* If a and b are specified in the fixed NB exlclusion list of the    */
        /* psf, set LJ and electrostatic interactions to zero.                */
        /* ================================================================== */
        for(int n=0; n<exNBN[k]; n++){
          if ((a == exNB[k][n].a && b== exNB[k][n].b) || (a == exNB[k][n].b && b== exNB[k][n].a)) {
            eps = 0.0;
            sig = 0.0;
            qq  = 0.0;
            #ifdef EWALD
            ljset[k].pot[a][b].exNB_flag = 1;
            #endif
            #ifdef SASA
            #ifndef BGO
            ljset[k].pot[a][b].pab = 0.0;
            #endif
            #ifdef BGO
            ljset[k].pot[a][b].pab = 0.58720;
            #endif
            #endif
          }
        }
        
	qq *= (ELEC * ELEQ * ELEQ * NA * 0.001);

	/* ================================================================== */
	/* Save the LJ interaction parameters to the ljset structure to be    */
	/* used in simulation.                                                */
	/* ================================================================== */
        ljset[k].pot[a][b].eps = eps;
        ljset[k].pot[a][b].sig = sig;
        ljset[k].pot[a][b].qq  = qq;
        //ljset[k].pot[a][b].id  = id;
      }
    }

   }//for k
#ifdef MPI
// printf("rank %d: ljset[0].pot[10][25].eps is %lf\n", mpi.my_rank, ljset[0].pot[10][25].eps); 
#endif
#endif //#ifndef STYPE

#ifdef STYPE
  /* ================================================================== */
  /* Loop around number of site types interactions and apply the mixing */
  /* rules.                                                             */
  /* ================================================================== */
  for(int a=0; a<NUM_ATOM; a++) {
    for(int b=0; b<NUM_ATOM; b++) {
      ljset[atom_type[a].atomid][atom_type[b].atomid].eps = 4.184 * sqrt(epstmp[0][atom_type[a].atomid] * epstmp[0][atom_type[b].atomid]);
      ljset[atom_type[a].atomid][atom_type[b].atomid].sig = 0.5 * (2.0 * sigtmp[0][atom_type[a].atomid] + 2.0 * sigtmp[0][atom_type[b].atomid]);//the 2.0 comes from the fact that CHARMM lists sig/2 
      ljset[atom_type[a].atomid][atom_type[b].atomid].id = 0; //0 means there are no special nbfix, native contacts,      etc.
      #ifdef S14    
      ljset14[atom_type[a].atomid][atom_type[b].atomid].eps = 4.184 * sqrt(epstmp[1][atom_type[a].atomid] * epstmp[1][atom_type[b].atomid]);
      ljset14[atom_type[a].atomid][atom_type[b].atomid].sig = 0.5 * (2.0 * sigtmp[1][atom_type[a].atomid] + 2.0 * sigtmp[1][atom_type[b].atomid]);//the 2.0 comes from the fact that CHARMM lists sig/2 
      #endif
    }
  }

  /* ================================================================== */
  /* If a and b are part of a NBFIX, use the specified interaction      */
  /* parameters.                                                        */
  /* ================================================================== */
  for (int n=0; n<NUM_NBFIX; n++){
    ljset[nbfix_prop[n].type1][nbfix_prop[n].type2].eps = nbfix_prop[n].eps*-4.184;
    ljset[nbfix_prop[n].type2][nbfix_prop[n].type1].eps = nbfix_prop[n].eps*-4.184;    
    ljset[nbfix_prop[n].type1][nbfix_prop[n].type2].sig = nbfix_prop[n].sig; 
    ljset[nbfix_prop[n].type2][nbfix_prop[n].type1].sig = nbfix_prop[n].sig; 
    //ljset[nbfix_prop[n].type1][nbfix_prop[n].type2].id = 10; 
    //ljset[nbfix_prop[n].type2][nbfix_prop[n].type1].id = 10; 
  }
	
  /* ================================================================== */
  /*                                                                    */
  /* Set up hash table of excluded interactions.                        */
  /*                                                                    */
  /* ================================================================== */
 int maxn=0;
 hashsize=1000;
  for(int k=0; k<sim.NB; k++){
	table[k] = (struct hashlist*) calloc (NPRIME+1, sizeof(struct hashlist));
		if(table[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for table\n"); exit(1);}

		for(int i=0; i<(NPRIME+1); i++){
			table[k][i].a = (int*) calloc(hashsize, sizeof(int));
			table[k][i].b = (int*) calloc(hashsize, sizeof(int));
			table[k][i].p = (int*) calloc(hashsize, sizeof(int));
			table[k][i].n = 0;
		}//for i
	
		for(int i=0; i<xintN[k]; i++){
			int a = xint[k][i].a;
			int b = xint[k][i].b;
			unsigned long int index = (a*b) % NPRIME;
			if(table[k][index].n == 0){
				table[k][index].a[0] = a;
				table[k][index].b[0] = b;
				table[k][index].p[0] = i;
			}
			else{
				int qa = table[k][index].n;
				table[k][index].a[qa] = a;
				table[k][index].b[qa] = b;
				table[k][index].p[qa]  = i;
			}
			table[k][index].n++;
			if( table[k][index].n>maxn)
      {
        maxn = table[k][index].n;
        if(maxn > hashsize){
          fprintf(stdout,"The hash table isn't big enough. Increase 'hashsize' in setlj.cpp\n");
        }
      }
		}//for i
  }//for k





#endif

}
