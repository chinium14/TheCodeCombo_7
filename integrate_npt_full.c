/* ======================================================================== */
/* integrate_npt_full.cpp                                                   */
/*                                                                          */
/*    This subroutine is one of several that integrate the equations of     */
/* motion.  The subroutine called is controled by sim.ID2 which is          */
/* specified in simul.input.  This subroutine does not do multiple time     */
/* step                                                                     */
/* See: Tuckerman, M., et. al., "Reversible multiple time scale molecular   */
/* dynamics," J. Chem. Phys., 97 (3) 1990-2001 1 August 1992.               */
/*                                                                          */
/* Passed Parameters:                                                       */
/*     ibox:  The box number                                                */
/*                                                                          */
/* ======================================================================== */
#ifdef PR_NPT
#include "defines.h"

  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
void forces(int);
void pbc_npt_full(int);
void integrate_nhcp_full(double, int);
void eigen(double a[][3], int n, double d[], double v[][3]);
void matrixmul(double a1[][3], double a2[][3], double a[][3], int);
void boxinv(int);
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void integrate_npt_full (int ibox)
{

  /* ================================================================== */
  /*                                                                    */
  /* Variables needed for this subroutine.                              */
  /*                                                                    */
  /* ================================================================== */
#define E2 1.0 / 6.0
#define E4 E2 / 20.0
#define E6 E4 / 42.0
#define E8 E6 / 72.0

  int k = ibox;
  double vtemp[3][3], btemp[3][3], veig[3], veigv[3][3];
  double ubox[3][3], tbox[3][3];
  int h;  //counter

  double AA, AA2[3], ARG2, BB[3], POLY;

  double u1, u2, u3, uv1, uv2, uv3;

  



  /* ================================================================== */
  /*                                                                    */
  /* Assign the time steps needed for integration.                      */
  /*                                                                    */
  /* ================================================================== */
  double dt	=	sim.dt;
  double dt_h	=	0.5* sim.dt;

 
  /* ================================================================== */
  /*                                                                    */
  /* Update the particle velocites, box velocites, thermostat           */
  /* velocities, and thermstat positions to dt/2.                       */
  /*                                                                    */
  /* ================================================================== */
  integrate_nhcp_full(dt , k);

  /* ================================================================== */
  /*                                                                    */
  /* Apply the Modified Velocity Verlet algorithm to advance the        */
  /* momenta and postions of the particles and the box.                 */
  /*                                                                    */
  /* ================================================================== */

  /* ----------------------------------------------- */
  /* Update the particle velocities by dt/2.         */
  /* ----------------------------------------------- */
	for(int i=0; i<box[k].boxns; i++){
	    double imass=100.0/pott[k][i].mas;
	    vv[k][i].x = vv[k][i].x + dt_h * imass *  ff[k][i].x;
        vv[k][i].y = vv[k][i].y + dt_h * imass *  ff[k][i].y;
        vv[k][i].z = vv[k][i].z + dt_h * imass *  ff[k][i].z;
	}

  /* ----------------------------------------------- */
  /* Obtain coefficients for positon and box update. */
  /* ----------------------------------------------- */
    h=0;
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			vtemp[i][j] = va[k][h];
			btemp[i][j] = axes[k][h];
			h++;
		}
	}	  
	eigen(vtemp, 3, veig, veigv);
	for(int i=0; i<3; i++){
		AA   = exp(dt_h * veig[i]);
		AA2[i]  = AA * AA;
		ARG2 = (veig[i] * dt_h) * (veig[i] * dt_h);
		POLY = (((E8 * ARG2 + E6) * ARG2 + E4) * ARG2 + E2) * ARG2 + 1.0;
		BB[i]   = AA * POLY * dt;
	}
  /* ----------------------------------------------- */
  /* Update the particle positions by dt.            */
  /* ----------------------------------------------- */
	for(int i=0; i<box[k].boxns; i++){
  		u1 = atom[k][i].x * veigv[0][0] + atom[k][i].y * veigv[1][0] + atom[k][i].z * veigv[2][0];
		u2 = atom[k][i].x * veigv[0][1] + atom[k][i].y * veigv[1][1] + atom[k][i].z * veigv[2][1];
		u3 = atom[k][i].x * veigv[0][2] + atom[k][i].y * veigv[1][2] + atom[k][i].z * veigv[2][2];
		uv1 = veigv[0][0] * vv[k][i].x + veigv[1][0] * vv[k][i].y + veigv[2][0] * vv[k][i].z;
		uv2 = veigv[0][1] * vv[k][i].x + veigv[1][1] * vv[k][i].y + veigv[2][1] * vv[k][i].z; 
		uv3 = veigv[0][2] * vv[k][i].x + veigv[1][2] * vv[k][i].y + veigv[2][2] * vv[k][i].z;
		u1 = u1 * AA2[0] + uv1 * BB[0];
		u2 = u2 * AA2[1] + uv2 * BB[1];
		u3 = u3 * AA2[2] + uv3 * BB[2];
		double ddx = u1 * veigv[0][0] + u2 * veigv[0][1] + u3 * veigv[0][2];
		double ddy = u1 * veigv[1][0] + u2 * veigv[1][1] + u3 * veigv[1][2];
		double ddz = u1 * veigv[2][0] + u2 * veigv[2][1] + u3 * veigv[2][2];
		atom[k][i].x = ddx;
		atom[k][i].y = ddy;
		atom[k][i].z = ddz;
		atnopbc[k][i].x = ddx;
		atnopbc[k][i].y = ddy;
		atnopbc[k][i].z = ddz;
	}

	


  /* ----------------------------------------------- */
  /* Update the box axes to dt.                      */
  /* ----------------------------------------------- */
	matrixmul(veigv, btemp, ubox, 2);

	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			ubox[i][j] = AA2[i] * ubox[i][j];
		}
	}

	matrixmul(veigv, ubox, tbox, 1);

	h=0;
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			axes[k][h] = tbox[i][j];
			h++;
		}
	}

	//force a symmetric box matrix
	axes[k][3] = axes[k][1];
	axes[k][6] = axes[k][2];
	axes[k][7] = axes[k][5];
	
	boxinv(k);

  /* ----------------------------------------------- */
  /* Apply periodic boundary conditions.             */
  /* ----------------------------------------------- */
	pbc_npt_full(k);

  
  /* ----------------------------------------------- */
  /* Update the forces.                              */
  /* ----------------------------------------------- */
	forces(k);

  /* ----------------------------------------------- */
  /* Update the velocities to dt/2.                  */
  /* ----------------------------------------------- */
       /* *********************************************** */
       /* Zero out the vvab accumulator.                  */
       /* *********************************************** */
	for(int i=0; i<9; i++) vvab[k][i] = 0.0;

	for (int i=0; i< box[k].boxns; i++){
	    double imass=100.0/pott[k][i].mas;
		vv[k][i].x = vv[k][i].x + dt_h * imass*  ff[k][i].x;
		vv[k][i].y = vv[k][i].y + dt_h * imass*  ff[k][i].y;
		vv[k][i].z = vv[k][i].z + dt_h * imass*  ff[k][i].z;
		uu[k][i].x = vv[k][i].x;
		uu[k][i].y = vv[k][i].y;
		uu[k][i].z = vv[k][i].z;
		vvab[k][0] = vvab[k][0] + vv[k][i].x * vv[k][i].x / imass;
		vvab[k][1] = vvab[k][1] + vv[k][i].x * vv[k][i].y / imass;
		vvab[k][2] = vvab[k][2] + vv[k][i].x * vv[k][i].z / imass;
		vvab[k][4] = vvab[k][4] + vv[k][i].y * vv[k][i].y / imass;
		vvab[k][5] = vvab[k][5] + vv[k][i].y * vv[k][i].z / imass;
		vvab[k][8] = vvab[k][8] + vv[k][i].z * vv[k][i].z / imass;
	}
  /* ----------------------------------------------- */
  /* Convert gm/mol to kg.  Also, imass is           */
  /* mulitiplied by 100 above so this is also        */
  /* cancelled here for vvab[k][i]                   */
  /* ----------------------------------------------- */

	    vvab[k][0] = vvab[k][0] / 10.0 / NA;
		vvab[k][1] = vvab[k][1] / 10.0 / NA;
		vvab[k][2] = vvab[k][2] / 10.0 / NA;
		vvab[k][4] = vvab[k][4] / 10.0 / NA;
		vvab[k][5] = vvab[k][5] / 10.0 / NA;
		vvab[k][8] = vvab[k][8] / 10.0 / NA;
		
	
	//to have a symmetric vvab matrix;
	vvab[k][3] = vvab[k][1];
	vvab[k][6] = vvab[k][2];
	vvab[k][7] = vvab[k][5];

  /* ================================================================== */
  /*                                                                    */
  /* Update the particle velocites, box velocites, thermostat           */
  /* velocities, and thermstat positions to dt/2.                       */
  /*                                                                    */
  /* ================================================================== */
  integrate_nhcp_full(dt , k);

}
#endif
