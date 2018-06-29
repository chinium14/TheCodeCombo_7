#ifdef MMDOS
#include "defines.h"


void vinit_hyb			(int);
void svalues			(int);
void pbc_all			(int);
void kinet				(int);
void calcvalue			(int);
void nblist				(int);
void integrate_vverlet	(int, double);
void flat_histogram		(int);
double ran2				(void);




void mmdos_hmc(int k){

	int		old_bin;
	int		new_bin;
	double  E_initial;
	double	E_final;
	double	U_initial;
	double	U_final;
	double	K_initial;
	double	K_final;
	double	vscal[NBOXES];
	double	beta_old;
	double	beta_new;
	double  arg=0.0;

	double  d_of_f = 	 ( 3.0 * box[k].boxns -2.0)* 0.5;
	double  scal   = 1.0/( 3.0 * box[k].boxns *0.5 * RG * 0.0010);

	vscal[k] = 1.0;								// initial scale for velocity

	
	E_initial	= en[k].totals;	
	U_initial	= en[k].potens;
	K_initial	= E_initial - U_initial;
	old_bin		=  (int) ((E_initial - sim_dos[k].e_begin)/sim_dos[k].e_width);
	sim.T[k]	= K_initial * scal;
	sim.kT[k]	= RG * sim.T[k]*0.001;
	beta_old	= 1.0/sim.kT[k];			
	beta_new	= beta_old;
	
	vinit_hyb(k);
	kinet(k); // what if this kE is not exactly the K_initial then en.total will change too
//	calcvalue(k);	// whether to call this calclvalue and update K_initial and E_initial	
	
	for(int l =0; l< box[k].boxns; l++){
		atom_temp[k][l]		= atom[k][l];		/* Back up coordinates with pdb			*/
		atnopbc_temp[k][l]	= atnopbc[k][l];	/* Back up coordinates without pdb		*/
		ff_temp[k][l]		= ff[k][l];			/* Back up all the force  components	*/
	}

	#ifdef PRESSURE
	pvir_temp[k]	= pvir[k];					/* Back up all the virial components	*/
	#endif
	en_temp[k]		= en[k];					/* Back up all the energy components	*/
				

	
	for (int j=0;j < sim_hyb.cyc_hybmd; j++){
		integrate_vverlet(k, vscal[k]);
		#ifdef DLIST
			nl_check(0,box[k].boxns,k);
			if(nl_flag[k] == 1) nblist(k);
		#endif

	}// j loop ends here
  
  /* ----------------------------------------------- */
  /* Calculate the properties of the system.         */
  /* ----------------------------------------------- */
	kinet(k);    											
	calcvalue(k);
	E_final = en[k].totals;
	U_final	= en[k].potens;
	K_final	= en[k].kinet;

	if (E_final >= sim_dos[k].e_begin && E_final < sim_dos[k].e_end){
		new_bin	=  (int) ((E_final - sim_dos[k].e_begin)/sim_dos[k].e_width);
		arg =	dos_hist[k][old_bin].g_of_e - dos_hist[k][new_bin].g_of_e + 
					d_of_f * (log(K_final/K_initial)) - beta_old*(K_final-K_initial);
	}
				
  /* ----------------------------------------------- */
  /* If the moves is accepted then accumulate the	 */
  /* calculated values for block averages and error  */
  /* bars else first return to original values and   */
  /* then accumulate statistics						 */
  /* ----------------------------------------------- */

	if (E_final >= sim_dos[k].e_begin && E_final < sim_dos[k].e_end && arg >0){
			dos_hist[k][new_bin].h_of_e ++;
			dos_hist[k][new_bin].k_of_e +=K_final;
			#ifdef NLIST
			if(sim_dos[k].dos_acc %sim.nlist==0) nblist(k); // call Nblist every nblist accepted moves
			#endif
			sim_dos[k].dos_acc ++;
	}//accepted
	
	else if (E_final >= sim_dos[k].e_begin && E_final < sim_dos[k].e_end && exp(arg)> ran2()){
			dos_hist[k][new_bin].h_of_e ++;
			dos_hist[k][new_bin].k_of_e +=K_final;
			#ifdef NLIST
			if(sim_dos[k].dos_acc %sim.nlist==0) nblist(k); // call Nblist every nblist accepted moves
			#endif
			sim_dos[k].dos_acc ++;
	}//accepted
	
	else{
			for(int l =0; l< box[k].boxns; l++){
				atom[k][l]		= atom_temp[k][l];		/* Back up coordinates with pdb			*/
				atnopbc[k][l]	= atnopbc_temp[k][l];	/* Back up coordinates without pdb		*/
				ff[k][l]		= ff_temp[k][l];		/* Back to old  force  components		*/
			}
			#ifdef PRESSURE
			pvir[k]		= pvir_temp[k];					/* Back to old  virial components		*/
			#endif
			en[k]		= en_temp[k];					/* Back to old  energy components		*/
			dos_hist[k][old_bin].h_of_e ++;
			dos_hist[k][old_bin].k_of_e += K_initial;
	}//rejected

	calcvalue(k); // this calcvalue may lead to a system with not K_initial but asssigned KE in case move rejected
	svalues(k);
}
#endif


