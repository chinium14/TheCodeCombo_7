/* ======================================================================== */
/* calcvalue.cpp                                                            */
/*                                                                          */
/*		This subroutine calculates the instantaneous values for the         */
/* density, volume, pressure, and temperature and assigns these values to   */
/* the "box[k]." variables.  It also takes the instantaneous potential      */
/* energies calculated in cbend(), cbond(), cnonbond(), and ctorsion() and  */
/* the instantaneious kinetic energy calculated in kinet() to calculate     */
/* the total instantaneous energy.                                          */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"
#ifdef PR_NPT
void boxinv(int);
#endif

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void calcvalue (int ibox)
{
  int k = ibox;
  
  /* ================================================================== */
  /*                                                                    */
  /* Calculate the total potential energy of the system using the       */
  /* contributions calculated during the force subroutines. Note:       */
  /* en[k].potens is the shifted energy.                                */
  /*                                                                    */
  /* ================================================================== */
#ifdef GOLIK
  en[k].poten  = en[k].bond + en[k].tors + en[k].nbond;
  en[k].potens = en[k].bond + en[k].tors + en[k].nbonds;
#endif
#ifndef GOLIK
#ifndef SASA
  en[k].poten  = en[k].bond + en[k].bend + en[k].tors + en[k].nbond + en[k].impr;
  en[k].potens = en[k].bond + en[k].bend + en[k].tors + en[k].nbonds + en[k].impr;
#endif
#ifdef SASA
  en[k].poten  = en[k].bond + en[k].bend + en[k].tors + en[k].nbond+ en[k].esasa + en[k].impr;
  en[k].potens = en[k].bond + en[k].bend + en[k].tors + en[k].nbonds+ en[k].esasa + en[k].impr;
#endif
#ifndef NEUTRAL
  en[k].poten  += en[k].urey;
  en[k].potens += en[k].urey;
#endif
#endif //GOLIK
  /* ================================================================== */
  /*                                                                    */
  /* Calculate the total energy (kinetic + potential) of the system.    */
  /* If the SASA option is defined, this energy also includes the       */
  /* solute-solvent interactions by adding the en[k].esasa term.        */
  /*                                                                    */
  /* ================================================================== */
  en[k].total  = en[k].poten  + en[k].kinet;
  en[k].totals = en[k].potens + en[k].kinet;

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the volume of the box.                                   */
  /*                                                                    */
  /* ================================================================== */
#ifndef PR_NPT
  box[k].vol = (box[k].lx * box[k].ly * box[k].lz);
#endif

#ifdef PR_NPT
  boxinv(k);
#endif

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the density of the box.                                  */
  /*                                                                    */
  /* ================================================================== */
  box[k].dens = (AU * 1.0E30) * (box[k].weight / box[k].vol);

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the temperature of the box.                              */
  /*                                                                    */
  /* ================================================================== */
  box[k].temp = (2.0 * en[k].kinet) / (box[k].nfree * RG * 1.0E-3);

  /* ================================================================== */
  /*                                                                    */
  /* If the PRESSURE option is defined, calculate the pressure and      */
  /* pressure tensor of the box.                                        */
  /*                                                                    */
  /*	First, the diagonal of the virial tensor, contained in the      */
  /* 0, 2, and 5 positions of the pvir[k].XXXX[i] variables, is stored  */
  /* in pvir[k].total.  The sum of the entire virial tensor is also     */
  /* calculated and stored in pvir[k].sum.  The pressure is calculated  */
  /* using the standard viral expression using the average of the       */
  /* diagonal entries of the virial tensor and stored in box[k].press.  */
  /* The entire pressure tensor is also calcualted and stored in        */
  /* box[k].cpress[i].													*/
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE

  /* ----------------------------------------- */
  /* The long-range correction.                */
  /* ----------------------------------------- */
  double lrcv = pvir[k].lrc / box[k].vol;

  /* ----------------------------------------- */
  /* Zero out the trace accumulators.          */
  /* ----------------------------------------- */
  pvir[k].tnbond = 0.0;
  pvir[k].tbond  = 0.0;
  pvir[k].tbend  = 0.0;
#ifndef NEUTRAL
  pvir[k].turey  = 0.0;
#endif
  pvir[k].ttors  = 0.0;
  pvir[k].timpr  = 0.0;
#ifdef SASA
  pvir[k].tsasa  = 0.0;
#endif
#ifdef EWALD
  pvir[k].tewald = 0.0;
#endif
  /* ----------------------------------------- */
  /* Determine the total virial tensor from    */
  /* from the contributions calculated in the  */
  /* force calculations.                       */
  /* ----------------------------------------- */
  for(int i=0; i<6; i++) {
	  pvir[k].sum[i] = pvir[k].nbond[i]+ pvir[k].bond[i] + pvir[k].bend[i] + pvir[k].tors[i] 	+pvir[k].impr[i]; 
#ifdef SASA 
	pvir[k].sum[i] += pvir[k].sasa[i]; 
#endif
#ifndef NEUTRAL
	pvir[k].sum[i] += pvir[k].urey[i];
#endif
#ifdef EWALD
	pvir[k].sum[i] += pvir[k].ewald_real[i] + pvir[k].ewald_intra[i] + pvir[k].ewald_kspace[i];
#endif //EWALD

    if( (i==0) || (i==2) || (i==5) ) {
	  pvir[k].sum[i] +=  lrcv/3.0;
      pvir[k].tnbond += pvir[k].nbond[i];
      pvir[k].tbond  += pvir[k].bond[i];
      pvir[k].tbend  += pvir[k].bend[i];
#ifndef NEUTRAL
	  pvir[k].turey  += pvir[k].urey[i];
#endif
      pvir[k].ttors  += pvir[k].tors[i];
	  pvir[k].timpr  += pvir[k].impr[i];
#ifdef SASA
	  pvir[k].tsasa  += pvir[k].sasa[i];
#endif
#ifdef EWALD
	  pvir[k].tewald += pvir[k].ewald_real[i] + pvir[k].ewald_intra[i] + pvir[k].ewald_kspace[i];
#endif //EWALD
	}
  }

  /* ----------------------------------------- */
  /* Calculate the trace of the virial tensor. */
  /* ----------------------------------------- */
  pvir[k].total = pvir[k].tnbond + pvir[k].tbond + pvir[k].tbend + pvir[k].ttors +pvir[k].timpr+ lrcv;
#ifdef SASA 
		pvir[k].total += pvir[k].tsasa; 
#endif 
#ifndef NEUTRAL
		pvir[k].total += pvir[k].turey;
#endif
#ifdef EWALD
		pvir[k].total += pvir[k].tewald; 

#ifdef EWALD_DEBUG
		FILE *sout2;
		char name2[100];


#ifdef MPI
		sprintf(name2,"./OUTPUT/BOX%d/ewald_virial%d.output",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
		sprintf(name2,"./OUTPUT/BOX%d/ewald_virial%d.output",k,k);
#endif
		sout2= fopen(name2,"a");
		fprintf(sout2, "%lf	%lf	%lf %lf	\n", en[k].tewald/3.0, pvir[k].tewald/3.0, en[k].tewald/3.0-pvir[k].tewald/3.0, en[k].rewald-en[k].rewalds);
		fclose(sout2);
#endif //EWALD_DEBUG
#endif //EWALD


  /* ----------------------------------------- */
  /* Calculate the box pressure. (kPa)         */
  /* ----------------------------------------- */
  box[k].press = (2.0*en[k].kinet+pvir[k].total)/(3.0*NA)*(1.0E30/box[k].vol);

  /* ----------------------------------------- */
  /* Calculate the pressure tensor components. */
  /* ----------------------------------------- */
  for(int i=0; i<6; i++){ 
    box[k].cpress[i] = (2.0*en[k].ken[i] + pvir[k].sum[i]) / NA * (1.0E30/box[k].vol);
  }
#endif
#ifdef ELASTIC
	elas[k].Pij[0] = box[k].cpress[0];
	elas[k].Pij[1] = box[k].cpress[2];
	elas[k].Pij[2] = box[k].cpress[5];
	elas[k].Pij[3] = box[k].cpress[4];
	elas[k].Pij[4] = box[k].cpress[3];
	elas[k].Pij[5] = box[k].cpress[1];
#endif

  /* ================================================================== */
  /*                                                                    */
  /* If the simulation mode is selected to use a Nose-Hoover            */
  /* thermostat, the extended Hamiltonian(the conserved quantity) is    */
  /* calculated.(kJ/mol)                                                */
  /*                                                                    */
  /* ================================================================== */
#ifndef PR_NPT
  if(sim.ID2 == 4 || sim.ID2 == 5 || sim.ID2 == 8)
  {
	double term1 = 0.0;
	double term2 = 0.0;

	for(int i=0; i<nhc[k].M; i++)
		term1 += nhc[k].zeta[i].v*nhc[k].zeta[i].v/2.0*nhc[k].Q*NA*10.0;
	
	for(int i=1; i<nhc[k].M; i++)
		term2 += RG*1E-3*sim.T[k]*nhc[k].zeta[i].r;

	en[k].H = en[k].totals+(double)box[k].nfree*RG*1E-3*sim.T[k]*nhc[k].zeta[0].r+term1 +term2;
  }
#endif
  /* ================================================================== */
  /*                                                                    */
  /* If a Parinello Rahman NPT simulation is being used, the conserved  */
  /* quantity is calculated.(kJ/mol)                                    */
  /*                                                                    */
  /* ================================================================== */
#ifdef PR_NPT
	ke_thermostat[k] = 0.0;
	sum_ktzeta_i[k] = 0.0;

	for(int i=0; i<nhc[k].M; i++)
		ke_thermostat[k] += nhc[k].zeta[i].v*nhc[k].zeta[i].v/2.0*nhc[k].Q*NA*10.0;
	
	for(int i=1; i<nhc[k].M; i++)
		sum_ktzeta_i[k] += RG*1E-3*sim.T[k]*nhc[k].zeta[i].r;
	
	cnf9ktzeta1[k]=((double)box[k].nfree+6.0)*RG*1E-3*sim.T[k]*nhc[k].zeta[0].r;
	
	pr_pv[k]= sim.P[k][0]*box[k].vol * NA * 1.0E-30;
	
	cka[k]=keaxes[k]*10.0 * NA;

	
	en[k].H = en[k].poten  + en[k].kinet + cnf9ktzeta1[k] + ke_thermostat[k] + sum_ktzeta_i[k]
		+ pr_pv[k] + cka[k];

#endif
  
  #ifdef TWHAM
  if(en[k].potens < pe_min[k]) pe_min[k] = en[k].potens;
  if(en[k].potens > pe_max[k]) pe_max[k] = en[k].potens;
  #endif

}
