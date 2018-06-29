/* ======================================================================== */
/* pme_calc.cpp                                                             */
/*            Calculate energy and force for PME                            */
/* ======================================================================== */

#include "defines.h"

#ifdef SPME

void charge_grid (int);
void pme_force (int);
#ifdef MC
void pme_force_mc (int);
#endif
void fourn(double[],unsigned long[],int,int);
double erfc (double);

void pme_calc (int k, double *energy)
{

  int ORDER = sim.order;
  int GRIDX = sim.grid[0];
  int GRIDY = sim.grid[1];
  int GRIDZ = sim.grid[2];

  double blx = box[k].lx;
  double bly = box[k].ly;
  double blz = box[k].lz;
  double vol = blx*bly*blz;

  double kap = sim.kappa[k];
  double fac = PI*PI/kap/kap;

  double wxx = 0.0;
  double wyx = 0.0;
  double wyy = 0.0;
  double wzx = 0.0;
  double wzy = 0.0;
  double wzz = 0.0;

  double enrec = 0.0;

  int ngrid = GRIDX*GRIDY*GRIDZ;
  int ggyz = GRIDY*GRIDZ;

  // calculate charge grid
  charge_grid(k);

  // back FFT
  fftw_execute(qgrid_back);
  
  for(int ind=1; ind<ngrid; ind++) {

    int k1 = ind/ggyz;
    int knd = ind - k1*ggyz;
    int k2 = knd/GRIDZ;
    int k3 = knd - k2*GRIDZ;


    int mx = k1  - (k1 > GRIDX/2 ? GRIDX: 0);
    int my = k2  - (k2 > GRIDY/2 ? GRIDY: 0);
    int mz = k3  - (k3 > GRIDZ/2 ? GRIDZ: 0);
    
    double mhatx = mx/blx;
    double mhaty = my/bly;
    double mhatz = mz/blz;
    
    double mm = mhatx*mhatx + mhaty*mhaty + mhatz*mhatz;
    double bb = bsp_modx[k1+1]*bsp_mody[k2+1]*bsp_modz[k3+1];
    double den = PI*vol*mm*bb;
    double eterm = exp(-fac*mm)/den;
    double wterm = -2.0*(fac*mm + 1.0)/mm;
    int index = index3D(k1,k2,k3);
    double qr = qgrida[index][0];
    double qi = qgrida[index][1];
    double struc2 = qr*qr + qi*qi;
    
    //fprintf(stderr,"%d %d %d  %lf  %lf  %lf\n",k1+1,k2+1,k3+1,mm,bb,eterm);
    //fprintf(stderr,"%d  %d  %d  \t  %d\n",k1,k2,k3,index);
    
    enrec += eterm*struc2;
    
    wxx += eterm*struc2*(wterm*mhatx*mhatx + 1.0);
    wyx += eterm*struc2*(wterm*mhaty*mhatx);
    wyy += eterm*struc2*(wterm*mhaty*mhaty + 1.0);
    wzx += eterm*struc2*(wterm*mhatz*mhatx);
    wzy += eterm*struc2*(wterm*mhatz*mhaty);
    wzz += eterm*struc2*(wterm*mhatz*mhatz + 1.0);
    
    qgrida[index][0] *= eterm;
    qgrida[index][1] *= eterm;

  }

  double fan = 0.5*QFACTOR/sim.epsRF[k];

  // Call the subroutines to calucate the k-space and 
  // self energy contributions to the electrostatic energy 
  energy[0] = enrec * fan;               // k-space ewald energy
  energy[1] = en[k].sewald;              // self-correction energy
  energy[2] = energy[0] - energy[1];     // k-space - self-energy

#ifdef PRESSURE
  pvir[k].ewald_kspace[0] = wxx*fan;
  pvir[k].ewald_kspace[1] = wyx*fan;
  pvir[k].ewald_kspace[2] = wyy*fan;
  pvir[k].ewald_kspace[3] = wzx*fan;
  pvir[k].ewald_kspace[4] = wzy*fan;
  pvir[k].ewald_kspace[5] = wzz*fan;
#endif

  // forward FFT
  fftw_execute(qgrid_forward);

  // calculate forces
  pme_force(k);


}


void pme_force (int k)
{

  int ORDER = sim.order;
  int GRIDX = sim.grid[0];
  int GRIDY = sim.grid[1];
  int GRIDZ = sim.grid[2];

  double blx = box[k].lx;
  double bly = box[k].ly;
  double blz = box[k].lz;
  double bhx = blx*0.5;  
  double bhy = bly*0.5;  
  double bhz = blz*0.5;  

  double kappa = sim.kappa[k];
  double kpi   = kappa / sqrt(PI);

  int nsite = box[k].boxns;

  double sumx = 0.0;
  double sumy = 0.0;
  double sumz = 0.0;

  for(int i=0; i<nsite-1; i++) {

    double fxi = 0.0;
    double fyi = 0.0;
    double fzi = 0.0;

    int xxo = (int)(ss[i].x) - ORDER;
    int yyo = (int)(ss[i].y) - ORDER;
    int zzo = (int)(ss[i].z) - ORDER;

    int xo = xxo;
    for(int xi=1; xi<=ORDER; xi++) {
      xo++;
      //int x = xo + 1 + (xo < 0 ? GRIDX : 0);
      int x = xo + (xo < 0 ? GRIDX : 0);
      double mxi  = mm[i][xi].x;
      double dmxi = dm[i][xi].x;

      int yo = yyo;
      for(int yi=1; yi<=ORDER; yi++) {
	yo++;
	//int y = yo + 1 + (yo < 0 ? GRIDY : 0);
	int y = yo + (yo < 0 ? GRIDY : 0);
	double myi  = mm[i][yi].y;
	double dmyi = dm[i][yi].y;
	
	int zo = zzo;
	for(int zi=1; zi<=ORDER; zi++) {
	  zo++;
	  //int z = zo + 1 + (zo < 0 ? GRIDZ : 0);
	  int z = zo + (zo < 0 ? GRIDZ : 0);
	  double mzi  = mm[i][zi].z;
	  double dmzi = dm[i][zi].z;
	  
  	  int index = index3D(x,y,z);
	  double term = qgrida[index][0];
	  
	  fxi -= dmxi *  myi *  mzi * term;
	  fyi -=  mxi * dmyi *  mzi * term;
	  fzi -=  mxi *  myi * dmzi * term;
	  
	}
      }
    }
    
    // NEED TO CHECK UNITS ********************************
    double qq = pott[k][i].qq;
    double fac = qq * QFACTOR / sim.epsRF[k];
    fxi *= fac*GRIDX/blx;
    fyi *= fac*GRIDY/bly;
    fzi *= fac*GRIDZ/blz;
    
    ff[k][i].x += fxi;
    ff[k][i].y += fyi;
    ff[k][i].z += fzi;
    
    sumx += fxi;
    sumy += fyi;
    sumz += fzi;
  
}
    
  int i = nsite-1;

  double fxi = 0.0;
  double fyi = 0.0;
  double fzi = 0.0;
  
  int xxo = (int)(ss[i].x) - ORDER;
  int yyo = (int)(ss[i].y) - ORDER;
  int zzo = (int)(ss[i].z) - ORDER;
  
  int xo = xxo;
  for(int xi=1; xi<=ORDER; xi++) {
    xo++;
    //int x = xo + 1 + (xo < 0 ? GRIDX : 0);
    int x = xo + (xo < 0 ? GRIDX : 0);
    double mxi  = mm[i][xi].x;
    double dmxi = dm[i][xi].x;
    
    int yo = yyo;
    for(int yi=1; yi<=ORDER; yi++) {
      yo++;
      //int y = yo + 1 + (yo < 0 ? GRIDY : 0);
      int y = yo + (yo < 0 ? GRIDY : 0);
      double myi  = mm[i][yi].y;
      double dmyi = dm[i][yi].y;
      
      int zo = zzo;
      for(int zi=1; zi<=ORDER; zi++) {
	zo++;
	//int z = zo + 1 + (zo < 0 ? GRIDZ : 0);
	int z = zo + (zo < 0 ? GRIDZ : 0);
	double mzi  = mm[i][zi].z;
	double dmzi = dm[i][zi].z;
	
	int index = index3D(x,y,z);
	double term = qgrida[index][0];
	
	fxi -= dmxi *  myi *  mzi * term;
	fyi -=  mxi * dmyi *  mzi * term;
	fzi -=  mxi *  myi * dmzi * term;
	
      }
    }
  }
  
  double qq = pott[k][i].qq;
  double fac = qq * QFACTOR / sim.epsRF[k];
  fxi *= fac*GRIDX/blx;
  fyi *= fac*GRIDY/bly;
  fzi *= fac*GRIDZ/blz;
  
  ff[k][i].x += fxi;
  ff[k][i].y += fyi;
  ff[k][i].z += fzi;
  
  sumx += fxi;
  sumy += fyi;
  sumz += fzi;
  

  // remove net force on atoms generated from the spline approximations
  sumx /= (double)(nsite);
  sumy /= (double)(nsite);
  sumz /= (double)(nsite);

  for(int i=0; i<nsite; i++) {
    ff[k][i].x -= sumx;
    ff[k][i].y -= sumy;
    ff[k][i].z -= sumz;
  }

}

#ifdef MC
void pme_calc_mc (int k, double *energy)
{

  int ORDER = sim.order;
  int GRIDX = sim.grid[0];
  int GRIDY = sim.grid[1];
  int GRIDZ = sim.grid[2];

  double blx = box[k].lx;
  double bly = box[k].ly;
  double blz = box[k].lz;
  double vol = blx*bly*blz;

  double kap = sim.kappa[k];
  double fac = PI*PI/kap/kap;

  double wxx = 0.0;
  double wyx = 0.0;
  double wyy = 0.0;
  double wzx = 0.0;
  double wzy = 0.0;
  double wzz = 0.0;

  double enrec = 0.0;

  int ngrid = GRIDX*GRIDY*GRIDZ;
  int ggyz = GRIDY*GRIDZ;

  // calculate charge grid
  charge_grid(k);

  // back FFT
  fftw_execute(qgrid_back);
  
  for(int ind=1; ind<ngrid; ind++) {

    int k1 = ind/ggyz;
    int knd = ind - k1*ggyz;
    int k2 = knd/GRIDZ;
    int k3 = knd - k2*GRIDZ;


    int mx = k1  - (k1 > GRIDX/2 ? GRIDX: 0);
    int my = k2  - (k2 > GRIDY/2 ? GRIDY: 0);
    int mz = k3  - (k3 > GRIDZ/2 ? GRIDZ: 0);
    
    double mhatx = mx/blx;
    double mhaty = my/bly;
    double mhatz = mz/blz;
    
    double mm = mhatx*mhatx + mhaty*mhaty + mhatz*mhatz;
    double bb = bsp_modx[k1+1]*bsp_mody[k2+1]*bsp_modz[k3+1];
    double den = PI*vol*mm*bb;
    double eterm = exp(-fac*mm)/den;
    double wterm = -2.0*(fac*mm + 1.0)/mm;
    int index = index3D(k1,k2,k3);
    double qr = qgrida[index][0];
    double qi = qgrida[index][1];
    double struc2 = qr*qr + qi*qi;
    
    //fprintf(stderr,"%d %d %d  %lf  %lf  %lf\n",k1+1,k2+1,k3+1,mm,bb,eterm);
    //fprintf(stderr,"%d  %d  %d  \t  %d\n",k1,k2,k3,index);
    
    enrec += eterm*struc2;
    
    wxx += eterm*struc2*(wterm*mhatx*mhatx + 1.0);
    wyx += eterm*struc2*(wterm*mhaty*mhatx);
    wyy += eterm*struc2*(wterm*mhaty*mhaty + 1.0);
    wzx += eterm*struc2*(wterm*mhatz*mhatx);
    wzy += eterm*struc2*(wterm*mhatz*mhaty);
    wzz += eterm*struc2*(wterm*mhatz*mhatz + 1.0);
    
    qgrida[index][0] *= eterm;
    qgrida[index][1] *= eterm;

  }

  double fan = 0.5*QFACTOR/sim.epsRF[k];

  // Call the subroutines to calucate the k-space and 
  // self energy contributions to the electrostatic energy 
  energy[0] = enrec * fan;               // k-space ewald energy
  energy[1] = en[k].sewald;              // self-correction energy
  energy[2] = energy[0] - energy[1];     // k-space - self-energy

#ifdef PRESSURE
  pvir[k].ewald_kspace[0] = wxx*fan;
  pvir[k].ewald_kspace[1] = wyx*fan;
  pvir[k].ewald_kspace[2] = wyy*fan;
  pvir[k].ewald_kspace[3] = wzx*fan;
  pvir[k].ewald_kspace[4] = wzy*fan;
  pvir[k].ewald_kspace[5] = wzz*fan;
#endif

  // forward FFT
  fftw_execute(qgrid_forward);

  // calculate forces
  pme_force_mc(k);


}


void pme_force_mc (int k)
{

  int ORDER = sim.order;
  int GRIDX = sim.grid[0];
  int GRIDY = sim.grid[1];
  int GRIDZ = sim.grid[2];

  double blx = box[k].lx;
  double bly = box[k].ly;
  double blz = box[k].lz;
  double bhx = blx*0.5;  
  double bhy = bly*0.5;  
  double bhz = blz*0.5;  

  double kappa = sim.kappa[k];
  double kpi   = kappa / sqrt(PI);

  int nsite = box[k].boxns;

  double sumx = 0.0;
  double sumy = 0.0;
  double sumz = 0.0;

  for(int i=0; i<nsite-1; i++) {

    double fxi = 0.0;
    double fyi = 0.0;
    double fzi = 0.0;

    int xxo = (int)(ss[i].x) - ORDER;
    int yyo = (int)(ss[i].y) - ORDER;
    int zzo = (int)(ss[i].z) - ORDER;

    int xo = xxo;
    for(int xi=1; xi<=ORDER; xi++) {
      xo++;
      //int x = xo + 1 + (xo < 0 ? GRIDX : 0);
      int x = xo + (xo < 0 ? GRIDX : 0);
      double mxi  = mm[i][xi].x;
      double dmxi = dm[i][xi].x;

      int yo = yyo;
      for(int yi=1; yi<=ORDER; yi++) {
	yo++;
	//int y = yo + 1 + (yo < 0 ? GRIDY : 0);
	int y = yo + (yo < 0 ? GRIDY : 0);
	double myi  = mm[i][yi].y;
	double dmyi = dm[i][yi].y;
	
	int zo = zzo;
	for(int zi=1; zi<=ORDER; zi++) {
	  zo++;
	  //int z = zo + 1 + (zo < 0 ? GRIDZ : 0);
	  int z = zo + (zo < 0 ? GRIDZ : 0);
	  double mzi  = mm[i][zi].z;
	  double dmzi = dm[i][zi].z;
	  
  	  int index = index3D(x,y,z);
	  double term = qgrida[index][0];
	  
	  fxi -= dmxi *  myi *  mzi * term;
	  fyi -=  mxi * dmyi *  mzi * term;
	  fzi -=  mxi *  myi * dmzi * term;
	  
	}
      }
    }
    
    // NEED TO CHECK UNITS ********************************
    double qq = pott[k][i].qq;
    double fac = qq * QFACTOR / sim.epsRF[k];
    fxi *= fac*GRIDX/blx;
    fyi *= fac*GRIDY/bly;
    fzi *= fac*GRIDZ/blz;
    
    ffox[k][i] += fxi;
    ffoy[k][i] += fyi;
    ffoz[k][i] += fzi;
    
    sumx += fxi;
    sumy += fyi;
    sumz += fzi;
  
}
    
  int i = nsite-1;

  double fxi = 0.0;
  double fyi = 0.0;
  double fzi = 0.0;
  
  int xxo = (int)(ss[i].x) - ORDER;
  int yyo = (int)(ss[i].y) - ORDER;
  int zzo = (int)(ss[i].z) - ORDER;
  
  int xo = xxo;
  for(int xi=1; xi<=ORDER; xi++) {
    xo++;
    //int x = xo + 1 + (xo < 0 ? GRIDX : 0);
    int x = xo + (xo < 0 ? GRIDX : 0);
    double mxi  = mm[i][xi].x;
    double dmxi = dm[i][xi].x;
    
    int yo = yyo;
    for(int yi=1; yi<=ORDER; yi++) {
      yo++;
      //int y = yo + 1 + (yo < 0 ? GRIDY : 0);
      int y = yo + (yo < 0 ? GRIDY : 0);
      double myi  = mm[i][yi].y;
      double dmyi = dm[i][yi].y;
      
      int zo = zzo;
      for(int zi=1; zi<=ORDER; zi++) {
	zo++;
	//int z = zo + 1 + (zo < 0 ? GRIDZ : 0);
	int z = zo + (zo < 0 ? GRIDZ : 0);
	double mzi  = mm[i][zi].z;
	double dmzi = dm[i][zi].z;
	
	int index = index3D(x,y,z);
	double term = qgrida[index][0];
	
	fxi -= dmxi *  myi *  mzi * term;
	fyi -=  mxi * dmyi *  mzi * term;
	fzi -=  mxi *  myi * dmzi * term;
	
      }
    }
  }
  
  double qq = pott[k][i].qq;
  double fac = qq * QFACTOR / sim.epsRF[k];
  fxi *= fac*GRIDX/blx;
  fyi *= fac*GRIDY/bly;
  fzi *= fac*GRIDZ/blz;
  
  ffox[k][i] += fxi;
  ffoy[k][i] += fyi;
  ffoz[k][i] += fzi;
  
  sumx += fxi;
  sumy += fyi;
  sumz += fzi;
  

  // remove net force on atoms generated from the spline approximations
  sumx /= (double)(nsite);
  sumy /= (double)(nsite);
  sumz /= (double)(nsite);

  for(int i=0; i<nsite; i++) {
    ffox[k][i] -= sumx;
    ffoy[k][i] -= sumy;
    ffoz[k][i] -= sumz;
  }

}
#endif //MC

#endif
