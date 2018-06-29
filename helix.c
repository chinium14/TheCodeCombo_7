/* ======================================================================== */
/* helix.cpp																*/
/*                                                                          */
/*		This subroutine calculates the helical content of a peptide         */
/*	based on its current phi and psi values and then return it				*/
/*                                                                          */
/* ======================================================================== */
#include "defines.h"
void helix(int k){

	int flag1		= 0;
	int flag2		= 0;
	double phi		= 0.0;
	double psi		= 0.0;
	int count		= 0;
	int helres		= 0;
	double helicity = 0.0;
#ifndef BETAPEP
	/* ----------------------------------------	*/
	/* The following is the definition of helix	*/
	/* ----------------------------------------	*/
	double phi_l_lim =-100.0*PI/180.0; 
	double phi_u_lim = -40.0*PI/180.0;
	double psi_l_lim = -67.0*PI/180.0;
	double psi_u_lim =  -7.0*PI/180.0;

	/* ----------------------------------------	*/
	/* Now compute fraction helicity           	*/
	/* ----------------------------------------	*/

	for(int i=0; i<torsN[k]; i++) {
		if(tors[k][i].phitag ==1){
			phi = tors[k][i].phi;
			flag1=1;
		}
		else if (tors[k][i].psitag ==1){
			psi = tors[k][i].psi;
			flag2=1;
		}
		if (flag1*flag2 !=0){
			if (phi> phi_l_lim && phi < phi_u_lim && psi> psi_l_lim && psi< psi_u_lim){
				helres++;
			}
			flag1=0; flag2=0;
			count ++;							// This keeps count of total number of residues -1
		}
		helicity = (double)helres/(double)count;
	}
#endif//NOBETAPEP
#ifdef BETAPEP
/* ----------------------------------------	*/
	/* The following is the definition of helix	*/
	/* ----------------------------------------	*/
	double deg_2_rad = PI/180.0;
	double a_width	= 20.0*deg_2_rad;
	double b_width	= 39.0*deg_2_rad;
	double del_width = 19.0*deg_2_rad;
	double hel_phi_psi = 0.0;
	double hel_phi 	= 0.0;
	double hel_psi	= 0.0;
	double phi_o	= -130.0*deg_2_rad;
	double psi_o	= -140.0*deg_2_rad;
/*
	double phi_l_lim =-160.0*PI/180.0; 
	double phi_u_lim = -110.0*PI/180.0;
	double psi_l_lim = -160.0*PI/180.0;
	double psi_u_lim =  -110.0*PI/180.0;
*/
	/* ----------------------------------------	*/
	/* Now compute fraction helicity           	*/
	/* ----------------------------------------	*/

	for(int i=0; i<torsN[k]; i++) {

		if(tors[k][i].phitag ==1){
			phi = tors[k][i].phi;
			flag1=1;
		}
		else if (tors[k][i].psitag ==1){
			psi = tors[k][i].psi;
			flag2=1;
		}
		if (flag1*flag2 !=0){
			double del_psi = fabs(psi-psi_o);
			double del_phi = fabs(phi-phi_o);

			if(del_psi < a_width) hel_psi=1.0;
			else if (del_psi > b_width) hel_psi =0.0;
			else hel_psi = 1.0-(del_psi-a_width)/del_width; 

			if(del_phi < a_width) hel_phi=1.0;
			else if (del_phi > b_width) hel_phi =0.0;
			else hel_phi = 1.0-(del_phi-a_width)/del_width;

				hel_phi_psi=hel_phi_psi + hel_phi*hel_psi; 
			flag1=0; flag2=0;
			count ++;							// This keeps count of total number of residues -1
		}
		helicity = hel_phi_psi/(double)count;
	}
#endif //BETAPEP
	
	ordparam[k].hel = helicity;
//	return(helicity);
}

double tors_av(int k, int tag){
	int count =0; double av_dih=0.0;
	if (tag==0){
		for(int i=0; i<torsN[k]; i++) {
			if(tors[k][i].thetag ==1){
				av_dih += tors[k][i].theta;
				count++;
			}
		}
		if (count==0) av_dih=0.0;
		else av_dih	= av_dih/count;
	}
	else if (tag==1){
		for(int i=0; i<torsN[k]; i++) {
			if(tors[k][i].psitag ==1){
				av_dih += tors[k][i].psi;
				count++;
			}
		}
		if (count==0) av_dih=0.0;
		else av_dih	= av_dih/count;
	}
	else if (tag==2){
		for(int i=0; i<torsN[k]; i++) {
			if(tors[k][i].phitag ==1){
				av_dih += tors[k][i].phi;
				count++;
			}
		}
		if (count==0) av_dih=0.0;
		else av_dih	= av_dih/count;
	}
	
	return(av_dih);
}
