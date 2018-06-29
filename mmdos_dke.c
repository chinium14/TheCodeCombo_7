#ifdef MMDOS
#include "defines.h"


void svalues			(int);
void kinet				(int);
void calcvalue			(int);
double ran2				(void);





void mmdos_dke(int k){

	int		old_bin;
	int		new_bin;
	double  E_initial;
	double	E_final;
	double	U_initial;
	double	K_initial;
	double	K_final;
	double  d_of_f = 	 ( 3.0 * box[k].boxns -2.0)* 0.5;

	
	E_initial	= en[k].totals;	
	U_initial	= en[k].potens;
	K_initial	= E_initial - U_initial;
	old_bin		=  (int) ((E_initial - sim_dos[k].e_begin)/sim_dos[k].e_width);
	
	#ifdef PRESSURE
	pvir_temp[k]	= pvir[k];								/* Back up all the virial components	*/
	#endif
	en_temp[k]		= en[k];								/* Back up all the energy components	*/
				

	double rand = (ran2()-0.5)*2.0;
	K_final		= K_initial + rand* sim_dos[k].d_max;
	E_final		= U_initial + K_final;
	en[k].kinet = K_final;
		
    if (E_final >= sim_dos[k].e_begin && E_final < sim_dos[k].e_end && K_final >0){
		new_bin	=  (int) ((E_final - sim_dos[k].e_begin)/sim_dos[k].e_width);
		double arg =	dos_hist[k][old_bin].g_of_e - dos_hist[k][new_bin].g_of_e +
						d_of_f * (log(K_final/K_initial));
	
		if (arg >0){
			dos_hist[k][new_bin].h_of_e ++;
			dos_hist[k][new_bin].k_of_e +=K_final; // No need to call Nblist
			sim_dos[k].dke_acc ++;
		}//accepted
		
		else if (exp(arg)> ran2()){
			dos_hist[k][new_bin].h_of_e ++;
			dos_hist[k][new_bin].k_of_e +=K_final;
			sim_dos[k].dke_acc ++;
		}//accepted
		
		else{
			#ifdef PRESSURE
			pvir[k]		= pvir_temp[k];					/* Back to old   components		*/
			#endif
			en[k]		= en_temp[k];					
			dos_hist[k][old_bin].h_of_e ++;
			dos_hist[k][old_bin].k_of_e += K_initial;
		}//rejected
	}// if K_final >0
	else{
		#ifdef PRESSURE
		pvir[k]		= pvir_temp[k];						/* Back to old  components		*/
		#endif
		en[k]		= en_temp[k];					
		dos_hist[k][old_bin].h_of_e ++;
		dos_hist[k][old_bin].k_of_e += K_initial;
	}// if K_final <0 or E out of bounds

	calcvalue(k); // updates en totals etc.
	svalues(k);
}
#endif

