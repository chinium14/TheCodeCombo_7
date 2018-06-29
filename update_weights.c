#ifdef FX_EDOS
/* This programm implemets Troyers idea of increasing round trips */
#include "defines.h"

void update_weights(int k)
{   //printf("hi");
	double sum_f,sum_e,sum_fe,sum_e2;
	double numo,deno,conv;
	int n_neig=2;//no of points taken on each side for regression 	
    int lm,rm;
	for(int i=0;i<sim_dos[k].e_bins;i++){
		sum_f=0.0;
		sum_e=0.0;
		sum_fe=0.0;
		sum_e2=0.0;
        numo=0.0;
		deno=0.0;
		if(i<n_neig)  {
			lm=i;
			rm=2*n_neig-lm;
		}//end if
		else if(i>(sim_dos[k].e_bins-n_neig-1)) {
			rm=sim_dos[k].e_bins-i-1;           
            lm=2*n_neig-rm;
			
		}//end elseif

		else {
			rm=n_neig;
			lm=n_neig;
		}// end else
			for(int j=i-lm;j<i+rm+1;j++) {
                sum_fe=sum_fe+dos_hist[k][j].f_of_e*(sim_dos[k].e_begin+sim_dos[k].e_width*j);
				sum_e=sum_e+sim_dos[k].e_begin+sim_dos[k].e_width*j;
				sum_f=sum_f+dos_hist[k][j].f_of_e;
				sum_e2=sum_e2+(sim_dos[k].e_begin+sim_dos[k].e_width*j)*(sim_dos[k].e_begin+sim_dos[k].e_width*j);
			}//end for
			numo=sum_fe-(sum_f)*(sum_e)/(2.0*n_neig+1);
			deno=sum_e2-(sum_e)*(sum_e)/(2.0*n_neig+1);
			dos_hist[k][i].m_of_e=numo/deno;
		}// end for

	//updating weights
            conv=0; double sum=0.0; double w_of_e0; int nw_of_e0;
		    for (int i=0;i<sim_dos[k].e_bins;i++){
				dos_hist_old[k][i].w_of_e	=	dos_hist[k][i].w_of_e;
				dos_hist[k][i].w_of_e		=	dos_hist[k][i].w_of_e+0.5*(log(fabs(dos_hist[k][i].m_of_e))
												-log(dos_hist[k][i].nw_of_e));
				sum		=	sum + exp(dos_hist[k][i].w_of_e-dos_hist[k][0].w_of_e);	
			}//end of for loop
			w_of_e0		=	dos_hist[k][0].w_of_e;
			nw_of_e0	=	dos_hist[k][0].nw_of_e;	
			for (int i=0;i<sim_dos[k].e_bins;i++){
 				dos_hist[k][i].w_of_e	=	dos_hist[k][i].w_of_e - w_of_e0 - log (sum);
				/* COmpute g of e from the old weights as the histogram nw correpond to those weights	*/
				dos_hist[k][i].g_of_e	=	log(dos_hist[k][i].nw_of_e) - dos_hist_old[k][i].w_of_e -log(nw_of_e0); 
				conv=conv+(dos_hist[k][i].w_of_e-dos_hist_old[k][i].w_of_e)*(dos_hist[k][i].w_of_e-dos_hist_old[k][i].w_of_e);
				dos_hist[k][i].nw_of_e	=	0;
				dos_hist[k][i].np_of_e	=	0;
				dos_hist[k][i].nm_of_e	=	0;

			}
			conv=conv/sim_dos[k].e_bins;
			sim_dos[k].tol_w=conv;


		}

#endif		
			
		






