
#include "defines.h"

#ifdef NMA
void forces(int);
void nblist(int);
void pbc_all(int);
void calcvalue(int);
#ifdef DLIST
void nl_check(int, int, int);
#endif
void force_nma(double *pos_nma,double *for_nma)
{
	for (int i=0; i<box[0].boxns; i++){
		atom[0][i].x = pos_nma[i*3+1];
		atom[0][i].y = pos_nma[i*3+2];
		atom[0][i].z = pos_nma[i*3+3];
	}
	pbc_all(0);
#ifdef NLIST
#ifndef DLIST
	nblist(0);
#endif
#ifdef DLIST
	nl_check(0, box[k].boxns,0);
	if(nl_flag[k] == 1) nblist(0);
#endif
#endif

	forces(0);
	for (int i=0; i<box[0].boxns; i++){
		for_nma[i*3+1] = ff[0][i].x;
		for_nma[i*3+2] = ff[0][i].y;
		for_nma[i*3+3] = ff[0][i].z;
	}
}
double energy_nma(double *pos_nma)
{
	for (int i=0; i<box[0].boxns; i++){
		atom[0][i].x = pos_nma[i*3+1];
		atom[0][i].y = pos_nma[i*3+2];
		atom[0][i].z = pos_nma[i*3+3];
	}
	pbc_all(0);
#ifdef NLIST
#ifndef DLIST
	nblist(0);
#endif
#ifdef DLIST
	nl_check(0, box[k].boxns,0);
	if(nl_flag[k] == 1) nblist(0);
#endif
#endif

	forces(0);
	calcvalue(0);
	return en[0].potens;
}
#endif
