/* ======================================================================== */
/* outend.cpp                                                               */
/*                                                                          */
/*		This subroutine writes the simulation statistics to the output      */
/* file.  The averages and root mean square deviations for each property is */
/* calculated and written. It is called only in the program, at the end.    */
/* This file also contains a subroutine called "rms" that calculates the    */
/* root mean square deviation of each property.                             */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
double rms (int, double, double);

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void outend (unsigned long icycle)
{

  /* ================================================================== */
  /* Calculate and write the total simulation time and production time. */
  /* ================================================================== */
  double totime_eq = sim.dt * sim.cyc_eq;
  double totime_pr;

  if (sim.ID2==3 || sim.ID2==5) totime_pr =sim.dtlong*1.0*icycle; // for multiple time step
  else totime_pr = sim.dt * 1.0*(icycle);
  
  /* ================================================================== */
  /*                                                                    */
  /* Calculate the average and rms of each property.  Note: These       */
  /* will be the average of the block averages and the rms of the same  */
  /* block averages.  In other words, it effectively averages over      */
  /* many different simulations.                                        */
  /*                                                                    */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++) {

#ifdef MPI
	  sprintf(name,"./OUTPUT/BOX%d/simul%d.output",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	 sprintf(name,"./OUTPUT/BOX%d/simul%d.output",k,k);
#endif  
	
	int nb = resb[k].count;
	
        //simul_ptr[k] += sprintf(simul_ptr[k]," Results - Ensemble Average:\n\n");
	//simul_ptr[k] += sprintf(simul_ptr[k],"Total equilibration time:         %12.4f \n",totime_eq);
	//simul_ptr[k] += sprintf(simul_ptr[k],"Total production time:            %12.4f \n",totime_pr);
	
	simul_ptr[k] += sprintf(simul_ptr[k],"Running Averages and Errors from %d Block(s) (blockd):\n\n",nb);
	simul_ptr[k] += sprintf(simul_ptr[k],"Total equilibration time:         %12.4f \n",totime_eq);
	simul_ptr[k] += sprintf(simul_ptr[k],"Total production time:            %12.4f \n",totime_pr);
	

  /* ----------------------------------------------- */
  /* Temperature                                     */
  /* ----------------------------------------------- */    
	double atemp = resb[k].tempb / (1.0*nb);
	double btemp = rms(nb,resb[k].tempb,resb[k].tempc);    
	simul_ptr[k] += sprintf(simul_ptr[k],"%d Temperature_av:                   %12.4f   %12.4f \n",k,atemp,btemp);

  /* ----------------------------------------------- */
  /* Density                                         */
  /* ----------------------------------------------- */    
    double adens = resb[k].densb / (1.0*nb);
    double bdens = rms(nb,resb[k].densb,resb[k].densc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Density_av:                       %12.4f   %12.4f \n",k,adens,bdens);

  /* ----------------------------------------------- */
  /* Anisotropic Pressure and Pressure Tensor        */
  /* ----------------------------------------------- */    
#ifdef PRESSURE
    double apress = resb[k].pressb / (1.0*nb);
    double bpress = rms(nb,resb[k].pressb,resb[k].pressc);
	simul_ptr[k] += sprintf(simul_ptr[k],"%d Pressure_av (atomic):             %12.4f   %12.4f \n",k,apress,bpress);


    
    double acpress[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    double bcpress[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    for(int i=0; i<6; i++) {
      acpress[i] = resb[k].cpressb[i] / (1.0*nb);
      bcpress[i] = rms(nb,resb[k].cpressb[i],resb[k].cpressc[i]);
	}
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Press_av(xx):                     %12.4f   %12.4f \n",k,acpress[0],bcpress[0]);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Press_av(yy):                     %12.4f   %12.4f \n",k,acpress[2],bcpress[2]);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Press_av(zz):                     %12.4f   %12.4f \n",k,acpress[5],bcpress[5]);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Press_av(xy):                     %12.4f   %12.4f \n",k,acpress[1],bcpress[1]);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Press_av(xz):                     %12.4f   %12.4f \n",k,acpress[3],bcpress[3]);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Press_av(yz):                     %12.4f   %12.4f \n",k,acpress[4],bcpress[4]);
#endif


  /* ----------------------------------------------- */
  /* Total energy (unshifted)                        */
  /* ----------------------------------------------- */    
    double aetotal = resb[k].etotalb / (1.0*nb);
    double betotal = rms(nb,resb[k].etotalb,resb[k].etotalc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Total energy_av:                  %12.4f   %12.4f \n",k,aetotal,betotal);

  /* ----------------------------------------------- */
  /* Total energy (shifted)                          */
  /* ----------------------------------------------- */    
    double aetotals = resb[k].etotalsb / (1.0*nb);
    double betotals = rms(nb,resb[k].etotalsb,resb[k].etotalsc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Total energy_av (shifted):        %12.4f   %12.4f \n",k,aetotals,betotals);

  /* ----------------------------------------------- */
  /* Kinetic Energy                                  */
  /* ----------------------------------------------- */    
    double aetkin = resb[k].etkinb / (1.0*nb);
    double betkin = rms(nb,resb[k].etkinb,resb[k].etkinc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Total kinetic energy_av:          %12.4f   %12.4f \n",k,aetkin,betkin);

  /* ----------------------------------------------- */
  /* Potential Energy (unshifted)                    */
  /* ----------------------------------------------- */    
    double aetpot = resb[k].etpotb / (1.0*nb);
    double betpot = rms(nb,resb[k].etpotb,resb[k].etpotc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Total potential energy_av:        %12.4f   %12.4f \n",k,aetpot,betpot);
#ifdef REM
	rep_exc[k].sig	= betpot;
	rep_exc[k].E_m	= aetpot;
#endif

  /* ----------------------------------------------- */
  /* Potential Energy (shifted)                      */
  /* ----------------------------------------------- */    
    double aetpots = resb[k].etpotsb / (1.0*nb);
    double betpots = rms(nb,resb[k].etpotsb,resb[k].etpotsc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Total pot. energy_av (shifted):   %12.4f   %12.4f \n",k,aetpots,betpots);

  /* ----------------------------------------------- */
  /* Nonbonded Energy (unshifted)                    */
  /* ----------------------------------------------- */    
    double aenbond = resb[k].enbondb / (1.0*nb);
    double benbond = rms(nb,resb[k].enbondb,resb[k].enbondc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Nonbonded energy_av:              %12.4f   %12.4f \n",k,aenbond,benbond);

  /* ----------------------------------------------- */
  /* Nonbonded Energy (shifted)                      */
  /* ----------------------------------------------- */    
    double aenbonds = resb[k].enbondsb / (1.0*nb);
    double benbonds = rms(nb,resb[k].enbondsb,resb[k].enbondsc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Nonbonded energy_av (shifted):    %12.4f   %12.4f \n",k,aenbonds,benbonds);

  /* ----------------------------------------------- */
  /* Bond Energy                                     */
  /* ----------------------------------------------- */    
    double aebond = resb[k].ebondb / (1.0*nb);
    double bebond = rms(nb,resb[k].ebondb,resb[k].ebondc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Bond energy_av:                   %12.4f   %12.4f \n",k,aebond,bebond);

  /* ----------------------------------------------- */
  /* Bend Energy                                     */
  /* ----------------------------------------------- */    
    double aebend = resb[k].ebendb / (1.0*nb);
    double bebend = rms(nb,resb[k].ebendb,resb[k].ebendc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Bend energy_av:                   %12.4f   %12.4f \n",k,aebend,bebend);

  /* ----------------------------------------------- */
  /* Urey-Bradley 1-3 Energy                         */
  /* ----------------------------------------------- */    
#ifndef NEUTRAL
    double aeurey = resb[k].eureyb / (1.0*nb);
    double beurey = rms(nb,resb[k].eureyb,resb[k].eureyc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Urey-Bradley 1-3 energy_av:       %12.4f   %12.4f \n",k,aeurey,beurey);
#endif

  /* ----------------------------------------------- */
  /* Torsion Energy                                  */
  /* ----------------------------------------------- */    
    double aetors = resb[k].etorsb / (1.0*nb);
    double betors = rms(nb,resb[k].etorsb,resb[k].etorsc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Torsion energy_av:                %12.4f   %12.4f \n",k,aetors,betors);

  /* ----------------------------------------------- */
  /* Improper torsion Energy                         */
  /* ----------------------------------------------- */    
    double aeimpr = resb[k].eimprb / (1.0*nb);
    double beimprs = rms(nb,resb[k].eimprb,resb[k].eimprc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Improper Torsion energy_av:       %12.4f   %12.4f \n",k,aeimpr,beimprs);

  /* ----------------------------------------------- */
  /* Heat Capacity                                   */
  /* ----------------------------------------------- */    
    double aheat_cap = resb[k].Cvb / (1.0*nb);
    double bheat_cap = rms(nb,resb[k].Cvb,resb[k].Cvc);
    simul_ptr[k] += sprintf(simul_ptr[k],"%d Heat Capacity_av:                    %12.6f   %12.4f \n",k,aheat_cap,bheat_cap);

  /* ----------------------------------------------- */
  /* Extended Hamiltonian Energy (Nose-Hoover and    */
  /* Parinello Rahman only)                          */
  /* ----------------------------------------------- */    
	if(sim.ID == 6 || ((sim.ID2 == 4 || sim.ID2 == 5 || sim.ID2==8) && sim.ID ==1)) 
	{
		double aH = resb[k].Hb/(1.0*nb);
		double bH = rms(nb,resb[k].Hb,resb[k].Hc);
		simul_ptr[k] += sprintf(simul_ptr[k],"%d Nose-Hoover Hamiltonian_av:       %12.4f   %12.4f \n\n",k,aH,bH);
	}

  }//k loop

}



/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
double rms (int n, double s1, double s2)
{

  double error = 0.0;

  if(n <= 1) {
    error = 0.0;
    return(error);
  }

  double in = 1.0/n;
  double arg = in * s2 - (in * s1)*(in * s1);

  error = sqrt(arg)/sqrt((double)n-1);

  return(error);

}
