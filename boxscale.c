/* ======================================================================== */
/* boxscale.cpp                                                             */
/*            Scales coordinates and box for constant-pressure simulations  */
/*              Method of Berendsen et al. 1984                             */
/* ======================================================================== */

#include "defines.h"

double pscal (double,double,int,int);

void boxscale (int ibox)
{

  int k = ibox;

#ifdef DPCONTROL
  sim.P[k][0] += sim.dpress[k][0];
  sim.P[k][1] += sim.dpress[k][1];
  sim.P[k][2] += sim.dpress[k][2];
  double scalex = pscal(box[k].cpress[0],sim.P[k][0],0,k);
  double scaley = pscal(box[k].cpress[2],sim.P[k][1],1,k);
  double scalez = pscal(box[k].cpress[5],sim.P[k][2],2,k);
#else
  sim.P[k][0] += sim.dpress[k][0];
  double scalex = pscal(box[k].press,sim.P[k][0],0,k);
  double scaley = scalex;
  double scalez = scalex;
#endif

  for(int i=0; i<box[k].boxns; i++) {
    atom[k][i].x *= scalex;
    atom[k][i].y *= scaley;
    atom[k][i].z *= scalez;
    atnopbc[k][i].x *= scalex;
    atnopbc[k][i].y *= scaley;
    atnopbc[k][i].z *= scalez;
  }

  box[k].lx *= scalex;
  box[k].ly *= scaley;
  box[k].lz *= scalez;
  box[k].hx	= box[k].lx / 2.0;
  box[k].hy	= box[k].ly / 2.0;
  box[k].hz	= box[k].lz / 2.0;
  box[k].vol = box[k].lx * box[k].ly *box[k].lz;

	double bd0 = box[k].hx;
	double bd1 = box[k].hy;
	double bd2 = box[k].hz;
	double bd = bd0;
	if((bd0 <= bd1) && (bd0 <= bd2)) bd = bd0;
	if((bd1 <= bd0) && (bd1 <= bd2)) bd = bd1;
	if((bd2 <= bd0) && (bd2 <= bd1)) bd = bd2;

	if (bd < sim.rc) {
		fprintf(stdout,"Error: The half boxlength (%lf) is smaller than the cutoff (%lf)\n", bd, sim.rc);
		exit(111);
	}
	 
#ifdef CELLIST
  // (CELLIST is the wrong keyword, if anyone ever sees this)
  //Cell list need to be updated for next energy calculations
#endif


}




double pscal (double press, double press0, int i, int k)
{
  double fr = 1.0/3.0;
#ifndef MTS	// I'm not sure MTS is the correct keyword, if anyone ever sees this
  double dt = sim.dt;
#endif
#ifdef MTS
  double dt = sim.dtlong;
#endif

  double tmp = (1.0 + (press - press0) * sim.Pcomp[k][i] * dt / sim.Ptau[k][i]);
  if(tmp < 0.0) {
    fprintf(stdout,"===> problem with pressure control %d (negative argument: %lf)\n",i,tmp);
    exit(1);
  }
  double scale = pow(tmp,fr);

  return(scale);
}
