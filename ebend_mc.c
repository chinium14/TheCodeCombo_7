#ifdef MC
#include "defines.h"


/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

double ebend_mc (int ibox)
{

  int k = ibox;
  /* ================================================================== */
  /*                                                                    */
  /* Zero out the energy accumulator                                    */
  /*                                                                    */
  /* ================================================================== */
  double ebend = 0.0;

  /* ================================================================== */
  /*                                                                    */
  /* Zero out the virial accumulators                                   */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  double wxx = 0.0;
  double wyx = 0.0;
  double wyy = 0.0;
  double wzx = 0.0;
  double wzy = 0.0;
  double wzz = 0.0;
#endif

  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the energy and force associated with each bend.  The     */
  /* loop is over the number of bends.                                  */
  /*                                                                    */
  /* ================================================================== */

  for(int i=0; i<bendN[k]; i++) {
  /* ----------------------------------------- */
  /* Assign the atoms involved in the bending. */
  /* ----------------------------------------- */
    int ia = bend[k][i].a;
    int ib = bend[k][i].b;
    int ic = bend[k][i].c;

  /* ----------------------------------------- */
  /* Calculate the distances between the atoms */
  /* involved in the bending.                  */
  /* ----------------------------------------- */
    double dax = atom[k][ia].x - atom[k][ib].x;
    double day = atom[k][ia].y - atom[k][ib].y;
    double daz = atom[k][ia].z - atom[k][ib].z;
    double dcx = atom[k][ic].x - atom[k][ib].x;
    double dcy = atom[k][ic].y - atom[k][ib].y;
    double dcz = atom[k][ic].z - atom[k][ib].z;

  /* ----------------------------------------- */
  /* Apply minimum image convention.           */
  /* ----------------------------------------- */
    if(dax >  hx) dax -= lx;
    else if(dax < -hx) dax += lx;
    if(day >  hy) day -= ly;
    else if(day < -hy) day += ly;
    if(daz >  hz) daz -= lz;
    else if(daz < -hz) daz += lz;


    if(dcx >  hx) dcx -= lx;
    else if(dcx < -hx) dcx += lx;
    if(dcy >  hy) dcy -= ly;
    else if(dcy < -hy) dcy += ly;
    if(dcz >  hz) dcz -= lz;
    else if(dcz < -hz) dcz += lz;

  /* ----------------------------------------- */
  /* Calculate the angle of the bend and then  */
  /* the energy using the CHARMM force field.  */
  /* ----------------------------------------- */
    double drai = 1.0 / sqrt(dax*dax + day*day + daz*daz);
    double erax = dax * drai; 
    double eray = day * drai; 
    double eraz = daz * drai; 
    double drci = 1.0 / sqrt(dcx*dcx + dcy*dcy + dcz*dcz);
    double ercx = dcx * drci; 
    double ercy = dcy * drci; 
    double ercz = dcz * drci; 

    double cosphi = erax*ercx + eray*ercy + eraz*ercz;
    double phi = acos(cosphi);
    double isinphi = 1.0 / sqrt(1.0 - cosphi*cosphi);

    double rr = phi - bend[k][i].angeq;
	double fr = -bend[k][i].krbend * rr*2.0;
    ebend += -fr * rr * 0.5;

    double irasin = drai * isinphi * fr;
    double ircsin = drci * isinphi * fr;

  /* ----------------------------------------- */
  /* Calculate the forces.                     */
  /* ----------------------------------------- */
    double fxa = (erax * cosphi - ercx) * irasin;
    double fya = (eray * cosphi - ercy) * irasin;
    double fza = (eraz * cosphi - ercz) * irasin;
    double fxc = (ercx * cosphi - erax) * ircsin;
    double fyc = (ercy * cosphi - eray) * ircsin;
    double fzc = (ercz * cosphi - eraz) * ircsin;

		ffox[k][ia] += fxa;
		ffoy[k][ia] += fya;
		ffoz[k][ia] += fza;
		ffox[k][ib] += -(fxa + fxc);
		ffoy[k][ib] += -(fya + fyc);
		ffoz[k][ib] += -(fza + fzc);
		ffox[k][ic] += fxc;
		ffoy[k][ic] += fyc;
		ffoz[k][ic] += fzc;


  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
    wxx += (fxa * dax + fxc * dcx);
    wyx += (fya * dax + fyc * dcx);
    wyy += (fya * day + fyc * dcy);
    wzx += (fza * dax + fzc * dcx);
    wzy += (fza * day + fzc * dcy);
    wzz += (fza * daz + fzc * dcz);
#endif

  }

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bending          */
  /* contribution to the virial tensor to pvir[k].bend[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pviro[k].bend[0] = wxx;
  pviro[k].bend[1] = wyx;
  pviro[k].bend[2] = wyy;
  pviro[k].bend[3] = wzx;
  pviro[k].bend[4] = wzy;
  pviro[k].bend[5] = wzz;
#endif

  return (ebend);

}
#endif//MC
