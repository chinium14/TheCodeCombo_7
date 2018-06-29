/* ======================================================================== */
/* pme_qgrid.cpp                                                            */
/*            Setup charge grid for PME                                     */
/* ======================================================================== */

#include "defines.h"

#ifdef SPME

void fillspline (double,double*,double*);


void charge_grid (int k)
{

  int ORDER = sim.order;
  int GRIDX = sim.grid[0];
  int GRIDY = sim.grid[1];
  int GRIDZ = sim.grid[2];

  double blx = box[k].lx;
  double bly = box[k].ly;
  double blz = box[k].lz;
  int nsite = box[k].boxns;

  for(int i=0; i<nsite; i++) {
    double qx = atom[k][i].x / blx; 
    double qy = atom[k][i].y / bly; 
    double qz = atom[k][i].z / blz;

    ss[i].x = GRIDX * (qx+0.5);
    ss[i].y = GRIDY * (qy+0.5);
    ss[i].z = GRIDZ * (qz+0.5);
  }

  double *mn, *dmn;
  mn = (double*) calloc(ORDER+1,sizeof(double));
  if(mn == NULL) {fprintf(stderr,"ERROR: memory, pmesetup: mn\n"); exit(1);}
  dmn = (double*) calloc(ORDER+1,sizeof(double));
  if(dmn == NULL) {fprintf(stderr,"ERROR: memory, pmesetup: dmn\n"); exit(1);}  


  for(int i=0; i<nsite; i++) {
    double ssx = ss[i].x;
    double wx = ssx - (int)(ssx);
    fillspline(wx,mn,dmn);

    for(int n=1; n<=ORDER; n++) {
      mm[i][n].x = mn[n];
      dm[i][n].x = dmn[n];
    }

    double ssy = ss[i].y;
    double wy = ssy - (int)(ssy);
    fillspline(wy,mn,dmn);

    for(int n=1; n<=ORDER; n++) {
      mm[i][n].y = mn[n];
      dm[i][n].y = dmn[n];
    }

    double ssz = ss[i].z;
    double wz = ssz - (int)(ssz);
    fillspline(wz,mn,dmn);

    for(int n=1; n<=ORDER; n++) {
      mm[i][n].z = mn[n];
      dm[i][n].z = dmn[n];
    }
  }

  free(mn); free(dmn);

  int ngrid = (GRIDX)*(GRIDY)*(GRIDZ);
  for(int i=0; i<ngrid; i++) { qgrida[i][0] = 0.0; qgrida[i][1] = 0.0; }

  for(int i=0; i<nsite; i++) {
    double qq = pott[k][i].qq;

    int xxo = (int)(ss[i].x) - ORDER;
    int yyo = (int)(ss[i].y) - ORDER;
    int zzo = (int)(ss[i].z) - ORDER;

    int xo = xxo;
    for(int xi=1; xi<=ORDER; xi++) {
      xo++;
      //int x = xo + 1 + (xo < 0 ? GRIDX : 0);
      int x = xo + (xo < 0 ? GRIDX : 0);
      double mxi = mm[i][xi].x;

      int yo = yyo;
      for(int yi=1; yi<=ORDER; yi++) {
	yo++;
	//int y = yo + 1 + (yo < 0 ? GRIDY : 0);
	int y = yo + (yo < 0 ? GRIDY : 0);
	double myi = mm[i][yi].y;
	
	int zo = zzo;
	for(int zi=1; zi<=ORDER; zi++) {
	  zo++;
	  //int z = zo + 1 + (zo < 0 ? GRIDZ : 0);
	  int z = zo + (zo < 0 ? GRIDZ : 0);
	  double mzi = mm[i][zi].z;

  	  int index = index3D(x,y,z);
	  if(index > ngrid)
	{
		fprintf(stdout, "OH NOES! index exceeds ngrid in pme_qgrid.c!\n");
		exit(1);
	}
	  //fprintf(stderr,"%d  %d  %d \t %d \n",x,y,z,index);
	  qgrida[index][0] += mxi*myi*mzi*qq;

	}
      }
    }
  }



}

#endif
