/* ======================================================================== */
/* mcmove_bkp.cpp                                                           */
/*     This subroutine backs up the coordinates, forces, energies, and      */
/* configurational temperature to the respective _temp variable so that a   */
/* Monte Carlo move can be performed.                                       */
/*                                                                          */
/* mcmove_rstr                                                              */
/*     This subroutine restores the coordinates, forces, energies, and      */
/* configurationtional temperature to the correct variables from the        */
/* corresponding _temp variables is a Monte Carlo move is rejected.         */
/*                                                                          */
/* Written by Thomas A. Knotts IV, 1 Sep 2005.                              */
/* ======================================================================== */
#include "defines.h"

#ifdef DLIST
void nl_check     (int,int,int);
#endif
void nblist_pivot (int,int);

void mcmove_bkp(int lb, int ub, int ibox){

  int k = ibox;

	for(int l=lb; l<ub; l++){
		atom_temp[k][l]	    = atom[k][l];				    
		atnopbc_temp[k][l]	= atnopbc[k][l];	
	}
	
  for(int l=0; l<box[k].boxns; l++) ff_temp[k][l] = ff[k][l];				      

  en_temp[k]		= en[k];  

  #ifdef PRESSURE
	pvir_temp[k]	= pvir[k];					      
  #endif
	
  #ifdef CONFIGT
  config_temp[k] = config[k];
  #endif

}


void mcmove_rstr(int lb, int ub, int ibox){

  int k = ibox;

  for(int l=lb; l<ub; l++){
	  atom[k][l]		= atom_temp[k][l];	 
	  atnopbc[k][l]	= atnopbc_temp[k][l];
	}

	for(int l=0; l< box[k].boxns; l++) ff[k][l]	= ff_temp[k][l];		

	en[k]		= en_temp[k];

	#ifdef PRESSURE
	pvir[k]		= pvir_temp[k];
	#endif

  #ifdef CONFIGT
  config[k] = config_temp[k];
  #endif

	#ifndef DLIST
	nblist_pivot(k,ub);
	#else
	nl_check(lb,ub,k);
	if(nl_flag[k] == 1) nblist_pivot(k,ub);
	#endif


}

