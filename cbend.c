/* ======================================================================== */
/* cbend.cpp                                                                */
/*                                                                          */
/*		This subroutine calculates the bending energies and forces.         */
/* It returns the value of the potential energy due to bending to be stored */
/* in en[k].ebend.  It is called from either forces() or force_short().     */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

#ifdef PR_NPT
void min_image_npt_full(int, double *dx, double *dy, double *dz);
#endif
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

double cbend (int ibox)
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
#ifndef PR_NPT
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
#endif
#ifdef PR_NPT
	min_image_npt_full(k, &dax, &day, &daz);
	min_image_npt_full(k, &dcx, &dcy, &dcz);
#endif

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


    ff[k][ia].x += fxa;
    ff[k][ia].y += fya;
    ff[k][ia].z += fza;
    ff[k][ib].x += -(fxa + fxc);
    ff[k][ib].y += -(fya + fyc);
    ff[k][ib].z += -(fza + fzc);
    ff[k][ic].x += fxc;
    ff[k][ic].y += fyc;
    ff[k][ic].z += fzc;

  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
    double sxx = (fxa * dax + fxc * dcx);
    double syx = (fya * dax + fyc * dcx);
    double syy = (fya * day + fyc * dcy);
    double szx = (fza * dax + fzc * dcx);
    double szy = (fza * day + fzc * dcy);
    double szz = (fza * daz + fzc * dcz);
    wxx += sxx;
    wyx += syx;
    wyy += syy;
    wzx += szx;
    wzy += szy;
    wzz += szz;
  /* ----------------------------------------- */
  /* If ZHOU is defined, accumulate the        */
  /* components of the stress tensor.          */
  /* ----------------------------------------- */
#ifdef ZHOU
	double third = 1.0 / 3;
        stresses[k][ia][0][0] += sxx * third;
        stresses[k][ia][1][0] += syx * third;
        stresses[k][ia][1][1] += syy * third;
        stresses[k][ia][2][0] += szx * third;
        stresses[k][ia][2][1] += szy * third;
        stresses[k][ia][2][2] += szz * third;

        stresses[k][ib][0][0] += sxx * third;
        stresses[k][ib][1][0] += syx * third;
        stresses[k][ib][1][1] += syy * third;
        stresses[k][ib][2][0] += szx * third;
        stresses[k][ib][2][1] += szy * third;
        stresses[k][ib][2][2] += szz * third;

        stresses[k][ic][0][0] += sxx * third;
        stresses[k][ic][1][0] += syx * third;
        stresses[k][ic][1][1] += syy * third;
        stresses[k][ic][2][0] += szx * third;
        stresses[k][ic][2][1] += szy * third;
        stresses[k][ic][2][2] += szz * third;


#endif  //ZHOU

#endif
  }

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bending          */
  /* contribution to the virial tensor to pvir[k].bend[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pvir[k].bend[0] = wxx;
  pvir[k].bend[1] = wyx;
  pvir[k].bend[2] = wyy;
  pvir[k].bend[3] = wzx;
  pvir[k].bend[4] = wzy;
  pvir[k].bend[5] = wzz;
#endif

  return (ebend);

}
