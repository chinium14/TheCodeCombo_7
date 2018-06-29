/* ======================================================================== */
/* mcmove_pivot.cpp                                                         */
/*     This move is a simple pivot.  It selects one of the solute molecules */
/* and then a site on that molecule.  It then pivot on one of the axiis.    */
/*                                                                          */
/* Written by Thomas Knotts, 1 Sep 2005.                                    */ 
/* ======================================================================== */

#include "defines.h"

void   svalues        (int);
void   pbc_all        (int);
void   calcvalue      (int);
void   nblist_pivot   (int,int);
double ran2           (void);
int	   ran_int	      (int,int);
double calc_dist      (int,int,int);     
void   mcmove_bkp     (int,int,int);
void   mcmove_rstr    (int,int,int);
#ifdef MC
void   forces_mc      (int,int,int,int);
#else
void   forces		      (int);
#endif

#ifdef DLIST
void   nl_check       (int,int,int);
#endif



void mcmove_pivot(int ibox)
{
	int k	=	ibox;

	/* ---------------------------------------------------- */
	/* First, pick a random solute molecule to move along   */
	/* the reaction coordinate.  Then determine the upper   */
	/* and lower bounds for the sites to move.              */
	/* ---------------------------------------------------- */
  int lb;
  int ub;
  int molecule;
  if(molec.Nsolute == 1){		
    lb = 0;
    ub = molec.lsite[0];
  }
  else{
    molecule = ran_int(0, molec.Nsolute);
    lb = molec.fsite[molecule];//first site of molecule
    ub = molec.lsite[molecule];//first site of next molecule
    if((ub-lb) < 5){
      lb = 0;
      ub = molec.lsite[0];
    }
  }
  int p;
  p = ran_int(lb,ub);		// pick the pivot site : pick Calpha or main_chain N atom only
  while (strcmp(atom[k][p].name,"CA")!=0 && strcmp(atom[k][p].name,"N")!=0 && 
    strcmp(atom[k][p].name,"P") !=0){
    p = ran_int(lb,ub);
	}

	int AXIS	=	ran_int(0,3);				// rotation about x, y or z axis
  #ifdef WALL
	int SIDE = 0;
  #else
	int SIDE	=	ran_int(0,2);				// amino or carboxy side with amino =0 and carboxy =1 side
  #endif	
	int ANGLE	=	ran_int(0,2);				// 0 is 90 1 is -90

  int start, end;

  if (SIDE==0){start = lb; end = p;}
  else {start = p+1; end = ub;}
  /* ====================================================== */
  /* Back up the necessary information so that the system   */
  /* can be restored if the move is rejected.               */
  /* ====================================================== */
  mcmove_bkp(start,end,k);
  int end_max = end;

  #ifdef DNA_GOLIK
  int start2, end2;
  int other_strand, lb2, ub2, index_shift;
  if(molecule == 0) other_strand = 1;
  else other_strand = 0;
  lb2 = molec.fsite[other_strand];//first site of molecule
  ub2 = molec.lsite[other_strand];//first site of next molecule
  index_shift = end - start;
  if(SIDE==0){
    start2 = ub2 - index_shift;
    end2 = ub2;
  }
  else{
    start2 = lb2;
    end2 = lb2 + index_shift;
  }
  /* ====================================================== */
  /* Back up the necessary information so that the system   */
  /* can be restored if the move is rejected.               */
  /* ====================================================== */
  mcmove_bkp(start2,end2,k);
  int start_min = start;
  if(end2 > end) end_max = end2;
  if(start2 < start) start_min = start2;
  #endif

  /* ====================================================== */
  /* Set the pre-move variables.                            */
  /* ====================================================== */
	double U_initial;
	double U_final;

  #ifdef MC
  #ifndef DNA_GOLIK
	forces_mc(start,end,0,k);
  #else
  forces_mc(start_min,end_max,0,k);
  #endif
  #endif
	calcvalue(k);
	
	double beta	= 1.0/sim.kT[k];			

  #ifdef MC		
	U_initial		= enmc[k].o_potens;
  #else		
	U_initial		= en[k].potens;
  #endif

  /* ====================================================== */
  /* Now propose the move.                                  */
  /* ====================================================== */
	double theta	=	ran2()*PI*mc_pivot.theta[k];
	double costheta	=	cos(theta);
	double sintheta	=	sin(theta);
	
	/*--------------------------------------*/
	/* INITIALIZATION OF PIVOT MATRICES		*/
	/*--------------------------------------*/

	int ii,jj, kk,ll;	
	for(ii=0; ii<3; ii++)
	{
		for(jj=0; jj<2; jj++)
		{
			for(kk=0; kk<3; kk++)
			{
				for(ll=0; ll<3; ll++)
					PIVOT[ii][jj][kk][ll]=0;
			}
		}
	}
	
	PIVOT[0][0][0][0]= 1.0;				// Rx (clockwise)
	PIVOT[0][0][1][1]= costheta;
	PIVOT[0][0][1][2]=-sintheta;
	PIVOT[0][0][2][1]= sintheta;
	PIVOT[0][0][2][2]= costheta;

	PIVOT[0][1][0][0]= 1.0;				// Rx (counter clockwise)
	PIVOT[0][1][1][1]= costheta;
	PIVOT[0][1][1][2]= sintheta;
	PIVOT[0][1][2][1]=-sintheta;
	PIVOT[0][1][2][2]= costheta;
	
	PIVOT[1][0][0][0]= costheta;
	PIVOT[1][0][0][2]= sintheta;		// Ry(clockwise)
	PIVOT[1][0][1][1]= 1.0;
	PIVOT[1][0][2][0]=-sintheta;
	PIVOT[1][0][2][2]= costheta;

	PIVOT[1][1][0][0]= costheta;
	PIVOT[1][1][0][2]=-sintheta;		// Ry(counter clockwise)
	PIVOT[1][1][1][1]= 1.0;
	PIVOT[1][1][2][0]= sintheta;
	PIVOT[1][1][2][2]= costheta;
	
	PIVOT[2][0][0][0]= costheta;
	PIVOT[2][0][0][1]= sintheta;		// Rz(clockwise)
	PIVOT[2][0][1][0]=-sintheta;
	PIVOT[2][0][1][1]= costheta;
	PIVOT[2][0][2][2]= 1.0;

	PIVOT[2][1][0][0]= costheta;
	PIVOT[2][1][0][1]=-sintheta;		// Rz(counter clockwise)
	PIVOT[2][1][1][0]= sintheta;
	PIVOT[2][1][1][1]= costheta;
	PIVOT[2][1][2][2]= 1.0;



  for (int i=start; i<end; i++){
    double pos_x	=	PIVOT[AXIS][ANGLE][0][0]*(atnopbc[k][i].x-atnopbc[k][p].x)+ 
                    PIVOT[AXIS][ANGLE][0][1]*(atnopbc[k][i].y-atnopbc[k][p].y)+
                    PIVOT[AXIS][ANGLE][0][2]*(atnopbc[k][i].z-atnopbc[k][p].z)+ atnopbc[k][p].x;
    double pos_y	=	PIVOT[AXIS][ANGLE][1][0]*(atnopbc[k][i].x-atnopbc[k][p].x)+
                    PIVOT[AXIS][ANGLE][1][1]*(atnopbc[k][i].y-atnopbc[k][p].y)+
                    PIVOT[AXIS][ANGLE][1][2]*(atnopbc[k][i].z-atnopbc[k][p].z)+ atnopbc[k][p].y;
    double pos_z	=	PIVOT[AXIS][ANGLE][2][0]*(atnopbc[k][i].x-atnopbc[k][p].x)+
                    PIVOT[AXIS][ANGLE][2][1]*(atnopbc[k][i].y-atnopbc[k][p].y)+
                    PIVOT[AXIS][ANGLE][2][2]*(atnopbc[k][i].z-atnopbc[k][p].z)+ atnopbc[k][p].z;
    atnopbc[k][i].x	=	pos_x;
    atnopbc[k][i].y	=	pos_y;
    atnopbc[k][i].z	=	pos_z;
    atom[k][i].x	=	atnopbc[k][i].x;
    atom[k][i].y	=	atnopbc[k][i].y;
    atom[k][i].z	=	atnopbc[k][i].z;
  }

  #ifdef DNA_GOLIK
  for (int i=start2; i<end2; i++){
    double pos_x	=	PIVOT[AXIS][ANGLE][0][0]*(atnopbc[k][i].x-atnopbc[k][p].x)+ 
                    PIVOT[AXIS][ANGLE][0][1]*(atnopbc[k][i].y-atnopbc[k][p].y)+
                    PIVOT[AXIS][ANGLE][0][2]*(atnopbc[k][i].z-atnopbc[k][p].z)+ atnopbc[k][p].x;
    double pos_y	=	PIVOT[AXIS][ANGLE][1][0]*(atnopbc[k][i].x-atnopbc[k][p].x)+
                    PIVOT[AXIS][ANGLE][1][1]*(atnopbc[k][i].y-atnopbc[k][p].y)+
                    PIVOT[AXIS][ANGLE][1][2]*(atnopbc[k][i].z-atnopbc[k][p].z)+ atnopbc[k][p].y;
    double pos_z	=	PIVOT[AXIS][ANGLE][2][0]*(atnopbc[k][i].x-atnopbc[k][p].x)+
                    PIVOT[AXIS][ANGLE][2][1]*(atnopbc[k][i].y-atnopbc[k][p].y)+
                    PIVOT[AXIS][ANGLE][2][2]*(atnopbc[k][i].z-atnopbc[k][p].z)+ atnopbc[k][p].z;
    atnopbc[k][i].x	=	pos_x;
    atnopbc[k][i].y	=	pos_y;
    atnopbc[k][i].z	=	pos_z;
    atom[k][i].x	=	atnopbc[k][i].x;
    atom[k][i].y	=	atnopbc[k][i].y;
    atom[k][i].z	=	atnopbc[k][i].z;
  }
  #endif
	/* ----------------------------------------------- */
	/* Apply periodic boundary conditions.             */
	/* ----------------------------------------------- */
	pbc_all(k);


  #ifndef DLIST
  nblist_pivot(k,end_max); // call Nblist
  #else
//  if (SIDE==0) nl_check(lb,p,k);
//  else nl_check(p+1,ub,k);
  nl_check(start,end,k);
  #ifdef DNA_GOLIK
  nl_check(start2, end2, k);
  #endif
  if(nl_flag[k] == 1) nblist_pivot(k,end_max);
  #endif


	/* ----------------------------------------------- */
	/* Update the forces.                              */
	/* ----------------------------------------------- */
  #ifdef MC
  #ifndef DNA_GOLIK
	forces_mc(start, end, 1,k);
  #else
  forces_mc(start_min, end_max, 1,k);
  #endif
  #else
	forces(k);
  #endif

	calcvalue(k);

  #ifdef MC
	U_final	= enmc[k].n_potens;
  #else
	U_final	= en[k].potens;
  #endif

  /* ====================================================== */
  /* Check to see if the move is accepted.                  */
  /* ====================================================== */
  double arg = (U_initial-U_final)*beta;
  if (arg>0) mc_pivot.pivot_acc[k]++; //accepted
	else if(exp(arg) > ran2()) mc_pivot.pivot_acc[k]++; //accepted
  else{
    mcmove_rstr(start,end,k); //rejected
    #ifdef DNA_GOLIK
    mcmove_rstr(start2,end2,k);
    #endif
  }
}
