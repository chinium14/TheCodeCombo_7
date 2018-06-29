/*------------------------------------------------------------- */
/*	This subroutine optimizes the flux associated with the swap	*/
/* moves by placing the boxes on most suitable Temperatures		*/
/* ------------------------------------------------------------ */

#ifdef REM
#include "defines.h"

	/* ----------------------------------------------------	*/
	/* The current sigma and E_m have already been computed	*/
	/* in the outend subroutine.							*/
	/* ----------------------------------------------------	*/
void rem_opt()
{
	double sum_sig=0.0;
	for (int k=0; k<sim.NB-1; k++){
		double del_E	= rep_exc[k+1].E_m - rep_exc[k].E_m;
		double del_sig	= rep_exc[k+1].sig - rep_exc[k].sig;
		rep_jnc[k].E_m /= rep_jnc[k].count; 
		rep_jnc[k].sig	= rep_exc[k].sig+(del_sig/del_E)*(rep_jnc[k].E_m- rep_exc[k].E_m);
		sum_sig			+= rep_jnc[k].sig;
		rep_jnc[k].flux	= (rep_jnc[k].count*sim_hyb.cyc_swap*2.0/sim.blockd);
		rep_opt.Jm		+= rep_jnc[k].flux;
		rep_jnc[k].rho	= rep_exc[k].rho;
	}
	rep_opt.Jm	/= (k-1);
	for (int k=0; k<sim.NB-1; k++){
		rep_opt.Je		+= fabs(rep_opt.Jm - rep_jnc[k].flux);
		rep_jnc[k].rho	= (1.0 + (sum_sig/rep_jnc[k].sig)*(rep_opt.Jm-rep_jnc[k].flux)/rep_opt.Je);
	}


}

void init_rem(){

	for (int k=0; k<sim.NB; k++){
		rep_exc[k].tag	= k;
		rep_exc[k].rho	= 1.0/sim.NB;
		rep_jnc[k].E_m	= 0.0;
		rep_jnc[k].count= 0;
		rep_opt.Jm		= 0.0;
		rep_opt.Je		= 0.0;
	}
}
#endif//REM