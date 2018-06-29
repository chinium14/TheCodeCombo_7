#include "defines.h"
void nblist		  (int);
double ran2		  (void);

int swap_box(int ibox)
{
	int k=ibox;
	double delta_beta	= 1.0/sim.kT[k]-1.0/sim.kT[k+1];
//	double delta_energy = en[k].totals-en[k+1].totals;
	double delta_energy = en[k].potens-en[k+1].potens;
	
	double accep_crit=exp(delta_beta*delta_energy);

		if(accep_crit > ran2())
		{
			for(int i =0; i< box[k].boxns; i++){
				atom_temp[k][i]		= atom[k][i];				/* Back up coordinates of box k			*/
				atnopbc_temp[k][i]	= atnopbc[k][i];			/* Back up coordinates of box k			*/
				atom[k][i]			= atom[k+1][i];				/* Transfer coordinates of k+1 to k		*/
				atnopbc[k][i]		= atnopbc[k+1][i];			/* Transfer coordinates of k+1 to k		*/
				atom[k+1][i]		= atom_temp[k][i];			/* Transfer coordinates of k to k+1		*/
				atnopbc[k+1][i]		= atnopbc_temp[k][i];		/* Transfer coordinates of k to k+1		*/
				ff_temp[k][i]		= ff[k][i];
				ff[k][i]			= ff[k+1][i];
				ff[k+1][i]			= ff_temp[k][i];
				vv_temp[k][i]		= vv[k][i];
				vv[k][i]			= vv[k+1][i];
				vv[k+1][i]			= vv_temp[k][i];
			}
			double v_scale_1		= sqrt(sim.T[k]/sim.T[k+1]);
			for(int i =0; i< box[k].boxns; i++){
				vv[k][i].x			*= v_scale_1;
				vv[k][i].y			*= v_scale_1;
				vv[k][i].z			*= v_scale_1;
				vv[k+1][i].x		/= v_scale_1;
				vv[k+1][i].y		/= v_scale_1;
				vv[k+1][i].z		/= v_scale_1;
				uu[k][i]			= vv[k][i];
				uu[k+1][i]			= vv[k+1][i];

			}


				
				/*======================================*/
				/* Tranfer data from box k to temp		*/
				/*======================================*/
#ifdef PRESSURE
				pvir_temp[k]		= pvir[k];					
#endif
				en_temp[k]			= en[k];					
				
				
				/*======================================*/
				/* Tranfer data from box k+1 to box k	*/
				/*======================================*/
#ifdef PRESSURE
				pvir[k]				= pvir[k+1];
#endif
				en[k]				= en[k+1];									
				#ifdef NLIST
					 nblist(k);									
				#endif
				/*======================================*/
				/* Tranfer data from temp to box k+1	*/
				/*======================================*/
#ifdef PRESSURE
				pvir[k+1]			= pvir_temp[k];	
#endif
				en[k+1]				= en_temp[k];						
				#ifdef NLIST
					 nblist(k+1);								
				#endif
					 
#ifdef REM
				int tag				= rep_exc[k].tag;
				rep_exc[k].tag		= rep_exc[k+1].tag;
				rep_exc[k+1].tag	= tag;
				rep_jnc[k].E_m		+=	((en[k].potens+en[k+1].potens)*0.5);
				rep_jnc[k].count	++;
#endif
			return 1;
		}
		else return 0;
		
}
