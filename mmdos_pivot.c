#ifdef MMDOS
/* ======================================================================== */
/* mmdos_pivot.cpp															*/
/* Warning:  As the code is written, pivot moves are only applied to the    */
/* first molecule (molecule 0) since the random sites are chosen from       */
/* integers between 0 and mol[0].Nsite.     								*/
/* ======================================================================== */
	
#include "defines.h"

void svalues	  (int);
void pbc_all	  (int);
void forces		  (int);
void calcvalue	  (int);
void nblist		  (int);
double ran2		  (void);
int	   ran_int	  (int,int);
void flat_histogram(int);


void mmdos_pivot(int ibox)
{
	int k			=	ibox;
	double theta	=	(ran2()+0.1)*PI*mc_pivot.theta[k];
	double costheta	=	cos(theta);
	double sintheta	=	sin(theta);
	
	// define all 30 enteries rest 24 are zeroes already set in init.cpp
	
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



	int p;
	p		=   ran_int(50,mol[0].Nsite-50);		// pick the pivot site : pick Calpha or main_chain N atom only
	while (strcmp(atom[k][p].name,"CA")!=0 && strcmp(atom[k][p].name,"N")!=0)p	=   ran_int(50,mol[0].Nsite-50);
	int AXIS	=	ran_int(0,3);				// rotation about x, y or z axis
	int SIDE	=	ran_int(0,2);				// amino or carboxy side with amino =0 and carboxy =1 side
	int ANGLE	=	ran_int(0,2);				// 0 is 90 1 is -90

	double U_initial;
	double K_initial;
	double E_initial;
	double U_final;
	double K_final;
	double E_final;
	double  d_of_f =	(3*box[k].boxns -2)*0.5;

	int old_bin;
	int new_bin;

//	calcvalue(k); //not needed
	for(int l =0; l< mol[0].Nsite; l++){
		atom_temp[k][l]	= atom[k][l];				/* Back up coordinates with pdb			*/
		atnopbc_temp[k][l]	= atnopbc[k][l];		/* Back up coordinates without pdb		*/
	}
	for(int l =0; l< box[k].boxns; l++){
		ff_temp[k][l]		= ff[k][l];				/* Back up all the force  components	*/
	}
#ifdef PRESSURE
		pvir_temp[k]	= pvir[k];					/* Back up all the virial components	*/
#endif
		en_temp[k]		= en[k];					/* Back up all the energy components	*/
		
		E_initial		= en[k].totals;
		K_initial		= en[k].kinet;				// same as E_initial - U_initial 
		U_initial		= en[k].potens;
		old_bin			=  (int) ((E_initial - sim_dos[k].e_begin)/sim_dos[k].e_width);

		
	/* Now propose a pivot move	*/


		if (SIDE==0){
			for (int i=0; i<p; i++){
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
		}
		else{	//if SIDE ==1
			for (int i=p+1; i<mol[0].Nsite; i++){
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
		}
		/* ----------------------------------------------- */
		/* Apply periodic boundary conditions.             */
		/* ----------------------------------------------- */
		pbc_all(k);

		/* ----------------------------------------------- */
		/* Update the forces.                              */
		/* ----------------------------------------------- */
		forces(k);
		calcvalue(k);
		
		U_final	= en[k].potens;
//		K_final	= E_initial - U_final;
//		en[k].kinet  = K_final;				// update kinetic energy in the energy structure
//		new_bin	= old_bin;
		
		E_final = E_initial+U_final-U_initial;
		K_final = K_initial;
		if (E_final >= sim_dos[k].e_begin && E_final < sim_dos[k].e_end){
			new_bin	=  (int) ((E_final - sim_dos[k].e_begin)/sim_dos[k].e_width);
			double arg = dos_hist[k][old_bin].g_of_e - dos_hist[k][new_bin].g_of_e;
			if (arg>0){
				dos_hist[k][new_bin].h_of_e ++;
				dos_hist[k][new_bin].k_of_e += K_final;
				#ifdef NLIST
				 nblist(k);					// call Nblist 
				#endif
				mc_pivot.pivot_acc[k] ++;
				calcvalue(k);				// updates en.totals etc
				svalues(k);
//				fprintf(stdout,"In box %d, pivot move with theta = %lf at site %d is accepted\n",k,(theta/PI),p);
			}//accepted
			else if (exp (arg) > ran2()){
				dos_hist[k][new_bin].h_of_e ++;
				dos_hist[k][new_bin].k_of_e += K_final;
				#ifdef NLIST
				nblist(k); // call Nblist 
				#endif
				mc_pivot.pivot_acc[k] ++;
				calcvalue(k);
				svalues(k);
//				fprintf(stdout,"In box %d, pivot move with theta = %lf at site %d is accepted\n",k,(theta/PI),p);
			}//accepted
			else{
				for(int l =0; l< mol[0].Nsite; l++){
					atom[k][l]		= atom_temp[k][l];	 /* Back up coordinates with pdb	*/
					atnopbc[k][l]	= atnopbc_temp[k][l];/* Back up coordinates without pdb	*/
					ff[k][l]		= ff_temp[k][l];	 /* Back to old  force  components	*/
				}
				#ifdef PRESSURE
				pvir[k]		= pvir_temp[k];				 /* Back to old  virial components	*/
				#endif
				en[k]		= en_temp[k];				 /* Back to old  energy components	*/
				dos_hist[k][old_bin].h_of_e ++;
				dos_hist[k][old_bin].k_of_e += K_initial;
				calcvalue(k);
				svalues(k);
			}//rejected due to arg
		}
		else{											/* that is if K_final <0			*/
			for(int l =0; l< mol[0].Nsite; l++){
				atom[k][l]		= atom_temp[k][l];		/* Back up coordinates with pdb		*/
				atnopbc[k][l]	= atnopbc_temp[k][l];	/* Back up coordinates without pdb	*/
			}
			for(int l =0; l< box[k].boxns; l++){
				ff[k][l]		= ff_temp[k][l];				/* Back up all the force  components	*/
			}
			
			#ifdef PRESSURE
			pvir[k]		= pvir_temp[k];					/* Back to old  virial components	*/
			#endif
			en[k]		= en_temp[k];					/* Back to old  energy components	*/
			dos_hist[k][old_bin].h_of_e ++;
			dos_hist[k][old_bin].k_of_e += K_initial;
			calcvalue(k);
			svalues(k);
		}//rejected coz K_final <=0
}
#endif

		
