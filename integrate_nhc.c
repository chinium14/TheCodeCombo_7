/* ======================================================================== */
/* integrate_nhc.cpp                                                        */
/*                                                                          */
/*      This subroutine integrates the thermostat postions, thermostat      */
/* velocities, and particle velocities using the Nose Hoover Chain method.  */
/* (This is the NHC propogator using the Louiville formalism.) It is only   */
/* done if a constant temperature simulation (sim.ID = 1 or 2) is           */
/* selected in simul.input with the Nose Hoover thermostat integration      */
/* option (sim.ID2 = 4 or 5).  At the end of the subroutine, half the       */
/* time step passed through the argument list is integrated.  The only      */
/* affect of this subroutine on the particles themselves is the scaling of  */
/* the particle velocities at the end of the subroutine.                    */
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

#include "defines.h"

  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
void kinet(int);

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void integrate_nhc(double time_step, int ibox)
{

  /* ================================================================== */
  /*                                                                    */
  /* Declare and assign the variables needed in the calculations.       */
  /*                                                                    */
  /* ================================================================== */	
	int h,i,j;
	double akin;
	double scale;

	int k = ibox;
	int M = nhc[k].M;
	int M1 = M-1;
	int M2 = M-2;
	double Q = nhc[k].Q;
	double Nf = (double)box[k].nfree;
	double NkT = Nf*KB*sim.T[k]*1E-4; //kg*ang^2/ps^2
	double kT = KB*sim.T[k]*1E-4; //kg*ang^2/ps^2
	double dt;

	double dt2;
	double dt4;
	double dt8;
	double AA;

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the kinetic energy.                                      */
  /*                                                                    */
  /* ================================================================== */
	kinet(k);
	akin = 2.0*en[k].kinet*1000.0*1E-4/NA;//conver kJ/mol to kg*ang^2/ps^2/molec (The factor of 2 comes from the fact that we need twice the kinetic energy)

  /* ================================================================== */
  /*                                                                    */
  /* Reset the scale accumulator.                                       */
  /*                                                                    */
  /* ================================================================== */
	scale = 1.0;
	
  /* ================================================================== */
  /*                                                                    */
  /* Start the multiple time step/fancy integration.                    */
  /*                                                                    */
  /* ================================================================== */
	for(h=0; h<nhc[k].Nc; h++){ //loop around number of substeps
		for(i=0; i<nhc[k].Nys; i++){ //loop for higher order terms
  
  /* ----------------------------------------------- */
  /* Times needed in algorithm.                      */
  /* ----------------------------------------------- */
			dt = time_step*nhc[k].w[i]/(double)nhc[k].Nc;//time step
			dt2 = dt/2.0;
			dt4 = dt/4.0;
			dt8 = dt/8.0;
  /* ----------------------------------------------- */
  /* Update the forces on the first thermostat.      */
  /* ----------------------------------------------- */
			nhc[k].zeta[0].G = (akin-NkT)/Q;
			//for(j=1; j<nhc[k].M; j++){ 
			//	nhc[k].zeta[j].G = 	(Q*nhc[k].zeta[j-1].v*nhc[k].zeta[j-1].v-kT)/Q;
			//}

  /* ----------------------------------------------- */
  /* Update the thermostat velocites to a quarter    */
  /* time step.                                      */
  /* ----------------------------------------------- */
			nhc[k].zeta[M1].v += nhc[k].zeta[M1].G*dt4;
			for(j=0; j<M1; j++){
				AA = exp(-dt8*nhc[k].zeta[M1-j].v);
				nhc[k].zeta[M2-j].v = nhc[k].zeta[M2-j].v*AA*AA + dt4*nhc[k].zeta[M2-j].G*AA;
			}

  /* ----------------------------------------------- */
  /* Accumulate the scaling factor.                  */
  /* ----------------------------------------------- */
			AA = exp(-dt2*nhc[k].zeta[0].v);
			scale = scale * AA;

  /* ----------------------------------------------- */
  /* Update the thermostat postions to half a time   */
  /* step.                                           */
  /* ----------------------------------------------- */
			for(j=0; j<M; j++){
				nhc[k].zeta[j].r += nhc[k].zeta[j].v*dt2;
			}
			
  /* ----------------------------------------------- */
  /* Update the forces on the first thermostat.      */
  /* ----------------------------------------------- */
			nhc[k].zeta[0].G = (scale*scale*akin-NkT)/Q;
			
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

		}//loop i
	}//loop h
		
  /* ----------------------------------------------- */
  /* Scale the particle velocities.                  */
  /* ----------------------------------------------- */
	for(i=0; i<box[k].boxns; i++){
		vv[k][i].x = vv[k][i].x*scale;
		vv[k][i].y = vv[k][i].y*scale;
		vv[k][i].z = vv[k][i].z*scale;
		uu[k][i].x = uu[k][i].x*scale;
		uu[k][i].y = uu[k][i].y*scale;
		uu[k][i].z = uu[k][i].z*scale;
	}

}

