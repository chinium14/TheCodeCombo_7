/* ======================================================================== */
/* integrate_nhcp_full.cpp                                                  */
/*                                                                          */
/*      This subroutine integrates the thermostat and box postions,         */
/* thermostat and box velocities and particle valocities                    */
/* velocities, and particle velocities using the Parinello Rahman barostat  */
/* thermostated to a Nose-Hoover chain thermostat.                          */
/*                                                                          */
/* See: Martyna et. al. "Explicit reversible integrators for extended       */
/* systems dynamics," Molecular Physics, Vol 87, No 5 1117-1157 (1996).     */
/*                                                                          */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                      time_step:	The whole time step of an interation.   */
/*                                  This value varies depending on if       */
/*                                  a multiple time step algorithm is used. */
/*                                                                          */
/* ======================================================================== */
#ifdef PR_NPT
#include "defines.h"

  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
void kinet(int);
void eigen(double a[][3], int n, double d[], double v[][3]);

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void integrate_nhcp_full(double time_step, int ibox)
{

  /* ================================================================== */
  /*                                                                    */
  /* Declare and assign the variables needed in the calculations.       */
  /*                                                                    */
  /* ================================================================== */	
	int h,i,j;
	double akin;

	int k = ibox;
	int M = nhc[k].M;
	int M1 = M-1;
	int M2 = M-2;
	double Q = nhc[k].Q;
	double W = Wg[k];
	double Nf = (double)box[k].nfree; 
	double N9kT = (Nf+6.0)*KB*sim.T[k]*1E-4; //kg*ang^2/ps^2
	double kT = KB*sim.T[k]*1E-4; //kg*ang^2/ps^2
	double dt;

	double P = sim.P[k][0]*1.0E-31;  //convert kPa to kg/(ang*ps^2)
	double V = box[k].vol;

	double dt2;
	double dt4;
	double dt8;
	double AA, BB;

	double trvava = 0.0;
	double trva;

	double sc1, sc2, sc3;
	double uv1, uv2, uv3;

  /* ----------------------------------------------- */
  /* Vectors and Matrices needed for box variables.  */
  /* ----------------------------------------------- */
	double vtemp[3][3], veig[3], veigv[3][3];
	
	
  /* ================================================================== */
  /*                                                                    */
  /* Calculate the kinetic energy.                                      */
  /*                                                                    */
  /* ================================================================== */
	kinet(k);
	akin = 2.0*en[k].kinet*1000.0*1E-4/NA;//conver kJ/mol to kg*ang^2/ps^2/molec (The factor of 2 comes from the fact that we need twice the kinetic energy)

  /* ================================================================== */
  /*                                                                    */
  /* Update the forces acting on the box and the thermostats.           */
  /*                                                                    */
  /* ================================================================== */
  /* ----------------------------------------------- */
  /* Times needed in algorithm.                      */
  /* ----------------------------------------------- */
	dt = time_step;
	dt2 = dt/2.0;
	dt4 = dt/4.0;
	dt8 = dt/8.0;
  /* ----------------------------------------------- */
  /* Forces on box axes.                             */
  /* ----------------------------------------------- */	
	for(i=0; i<9; i++){
		Ga[k][i] = (vvab[k][i] + frab[k][i] + delta[i]*(akin/Nf-P*V))/W;
		trvava = trvava + va[k][i] * va[k][i];	
	}

  /* ----------------------------------------------- */
  /* Update the forces on the first thermostat.      */
  /* ----------------------------------------------- */
	nhc[k].zeta[0].G = (akin + W * trvava - N9kT)/Q;	

  /* ----------------------------------------------- */
  /* Update the thermostat velocites by dt/4.        */
  /* ----------------------------------------------- */
	nhc[k].zeta[M1].v += nhc[k].zeta[M1].G*dt4;
	for(i=0; i<M1; i++){
		AA = exp(-dt8*nhc[k].zeta[M1-i].v);
		nhc[k].zeta[M2-i].v = nhc[k].zeta[M2-i].v*AA*AA 
							+ dt4*nhc[k].zeta[M2-i].G*AA;
	}

  /* ----------------------------------------------- */
  /* Update the velocity of the axes by dt/4.        */
  /* ----------------------------------------------- */
	BB = exp(-nhc[k].zeta[0].v *dt8);
	for(i=0; i<9; i++){
		va[k][i] = va[k][i] * BB;
		va[k][i] = va[k][i] + Ga[k][i] * dt4;
		va[k][i] = va[k][i] * BB;
	}
  /* ----------------------------------------------- */
  /* Update the thermostat postions to dt/2.         */
  /* ----------------------------------------------- */
	for(i=0; i<M; i++){
		nhc[k].zeta[i].r += nhc[k].zeta[i].v*dt2;
	}

  /* ----------------------------------------------- */
  /* Scale the velocities of the particles.          */
  /* ----------------------------------------------- */
	trva = (va[k][0] + va[k][4] + va[k][8]) / Nf;
	
	h=0;
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			vtemp[i][j] = va[k][h] + delta[h]*(trva + nhc[k].zeta[0].v);
			//vtemp[i][j] = va[k][h] + trva + nhc[k].zeta[0].v;
			h++;
		}
	}

	eigen(vtemp, 3, veig, veigv);
	sc1 = exp(-veig[0] * dt2);
	sc2 = exp(-veig[1] * dt2);
	sc3 = exp(-veig[2] * dt2);

	for(i=0; i<9; i++) vvab[k][i] = 0.0;
	
	for(i=0; i<box[k].boxns; i++){
		uv1 = vv[k][i].x * veigv[0][0] + vv[k][i].y * veigv[1][0] + vv[k][i].z * veigv[2][0];
		uv2 = vv[k][i].x * veigv[0][1] + vv[k][i].y * veigv[1][1] + vv[k][i].z * veigv[2][1];
		uv3 = vv[k][i].x * veigv[0][2] + vv[k][i].y * veigv[1][2] + vv[k][i].z * veigv[2][2];
		
		uv1 = uv1 * sc1;
		uv2 = uv2 * sc2;
		uv3 = uv3 * sc3;

		vv[k][i].x = uv1 * veigv[0][0] + uv2 * veigv[0][1] + uv3 * veigv[0][2];
		vv[k][i].y = uv1 * veigv[1][0] + uv2 * veigv[1][1] + uv3 * veigv[1][2];
		vv[k][i].z = uv1 * veigv[2][0] + uv2 * veigv[2][1] + uv3 * veigv[2][2];
		
		uu[k][i].x = vv[k][i].x;
		uu[k][i].y = vv[k][i].y;
		uu[k][i].z = vv[k][i].z;
	
		double mass = pott[k][i].mas;

		vvab[k][0] = vvab[k][0] + vv[k][i].x * vv[k][i].x * mass;
		vvab[k][1] = vvab[k][1] + vv[k][i].x * vv[k][i].y * mass;
		vvab[k][2] = vvab[k][2] + vv[k][i].x * vv[k][i].z * mass;
		vvab[k][4] = vvab[k][4] + vv[k][i].y * vv[k][i].y * mass;
		vvab[k][5] = vvab[k][5] + vv[k][i].y * vv[k][i].z * mass;
		vvab[k][8] = vvab[k][8] + vv[k][i].z * vv[k][i].z * mass;
	}
  /* ----------------------------------------------- */
  /* Convert gm/mol to kg.                           */
  /* ----------------------------------------------- */
	    vvab[k][0] = vvab[k][0] / 1000.0 / NA;
		vvab[k][1] = vvab[k][1] / 1000.0 / NA;
		vvab[k][2] = vvab[k][2] / 1000.0 / NA;
		vvab[k][4] = vvab[k][4] / 1000.0 / NA;
		vvab[k][5] = vvab[k][5] / 1000.0 / NA;
		vvab[k][8] = vvab[k][8] / 1000.0 / NA;



	//to have a symmetric vvab matrix;
	vvab[k][3] = vvab[k][1];
	vvab[k][6] = vvab[k][2];
	vvab[k][7] = vvab[k][5];

  /* ----------------------------------------------- */
  /* Since the velocites have been rescaled, update  */
  /* the kinetic energy.                             */
  /* ----------------------------------------------- */
	kinet(k);
	akin = 2.0*en[k].kinet*1000.0*1E-4/NA;//conver kJ/mol to kg*ang^2/ps^2/molec (The factor of 2 comes from the fact that we need twice the kinetic energy)

  /* ----------------------------------------------- */
  /* Update the force on the box axes.               */
  /* ----------------------------------------------- */	
//	for(i=0; i<9; i++)Ga[k][i] = (vvab[k][i] + frab[k][i] + delta[i]*(akin/Nf-P*V))/W;
		
  /* ----------------------------------------------- */
  /* Update the velocity of the axes by dt/4.        */
  /* ----------------------------------------------- */
	trvava = 0.0;
	//BB = exp(-nhc[k].zeta[0].v *dt8);
	for(i=0; i<9; i++){
		Ga[k][i] = (vvab[k][i] + frab[k][i] + delta[i]*(akin/Nf-P*V))/W;
		va[k][i] = va[k][i] * BB;
		va[k][i] = va[k][i] + Ga[k][i] * dt4;
		va[k][i] = va[k][i] * BB;
		trvava = trvava + va[k][i] * va[k][i];
	}

  /* ----------------------------------------------- */
  /* Update the forces on the first thermostat.      */
  /* ----------------------------------------------- */
	nhc[k].zeta[0].G = (akin + W * trvava - N9kT)/Q;

  /* ----------------------------------------------- */
  /* Update the thermostat velocites by dt/4.        */
  /* ----------------------------------------------- */
  /* ----------------------------------------------- */
  /* Update the forces on the other thermostats, and */
  /* and the velocities of all the thermostats to a  */
  /* quarter time step.                              */
  /* ----------------------------------------------- */
	for(j=0; j<M-1; j++){
		AA = exp(-dt8*nhc[k].zeta[j+1].v);
		nhc[k].zeta[j].v = nhc[k].zeta[j].v*AA*AA+dt4*nhc[k].zeta[j].G*AA;
		nhc[k].zeta[j+1].G = (Q*nhc[k].zeta[j].v*nhc[k].zeta[j].v-kT)/Q;
	}
	nhc[k].zeta[M-1].v += nhc[k].zeta[M-1].G*dt4;

  /* ----------------------------------------------- */
  /* Calculate the kinetic energy of the axes for    */
  /* the conserved quantity.                         */
  /* ----------------------------------------------- */
	keaxes[k] = 0.5 * W * trvava;

}//end subroutine

#endif


