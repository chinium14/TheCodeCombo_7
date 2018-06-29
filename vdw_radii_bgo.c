#ifdef SASA
#ifdef BGO
#include "defines.h"
/* ======================================================================== */
/* ======================================================================== */
//void set_go_params(int ibox);
//int is_bond(int,int,int);
void vdw_radii_bgo(void){
//  FILE * hi;
//  hi=fopen("vdw_radii_bgo.dat","w");
/*
{"ABC",                "GLY","ALA","SER","CYS","VAL","THR","ILE","PRO","MET","ASP",
                       "ASN","LEU","LYS","GLU","GLN","ARG","HIS","PHE","TYR","TRP",

 "CYT","GUA","ADE","THY","HOH","SOD","CLA","LIP","CHO","EAM","POT","TRE", "GOL","F00", "F01","F02","F03","F04","F05","F06","F07","F08","F09","F10","F11","MEO"}
*/

  scale_factor_b = 1.0; 
  double vdw[47]={0.0, 3.86518, 4.01256, 4.12788, 3.74826, 4.28142, 4.15639, 4.57110, 4.51078, 4.49251, 4.27733, 
                       4.28348, 4.55731, 4.16101, 4.29564, 4.35527, 4.62059, 4.32923, 4.69671, 4.51586, 5.05033, 0.0,
                       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double vdw_measured[47]={0.0, 4.5417, 5.5491, 6.2183, 6.3746, 7.0256, 6.8406, 7.7378, 6.6179, 8.4121, 7.3823,
                                7.3412, 8.0562, 9.0014, 8.4108, 8.4633, 10.1201, 8.6180, 9.1289, 9.7016, 9.9918, 0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
// Original p_list values by POPS
/*  double p_list[47]={0.0, 1.05203, 0.99042, 0.98230, 0.95123, 0.98227, 0.89542, 1.04589, 0.98890, 0.97725, 0.88912,
                          0.87896, 1.02160, 0.65241, 0.77849, 0.79978, 0.71552, 0.79982, 1.02315, 0.86185, 1.01999, 0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
*/
  double p_list[47]={0.0, 1.05203, 0.99042, 0.98230, 0.95123, 0.98227, 0.89542, 1.04589, 0.98890, 0.97725, 0.88912,
                          0.87896, 1.02160, 0.67000, 0.82000, 0.79978, 0.79000, 0.79982, 1.02315, 0.93000, 1.01999, 0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double sasa_ref_sc[47]={0.0, 0.0, 108.3, 114.9, 128.6, 147.1, 135.7, 164.7, 132.7, 164.7, 133.7, 
                             138.7, 164.6, 174.6, 150.8, 155.4, 186.0, 159.2, 174.6, 179.9, 198.7, 0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double sasa_ref_bb[47]={0.0, 85.0, 62.5, 60.9, 57.7, 53.8, 56.2, 50.3, 56.9, 50.3, 56.7,
                               55.6, 50.3, 48.3, 53.0, 52.1, 46.2, 51.4, 48.3, 47.3, 43.8, 0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double sasa_m_sc[47]={0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double sasa_m_bb = 0.0; 
  double sasa_b_sc[47] ={0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double sasa_b_bb = 0.0;

  if (solvent_type == 1){
        double sasa_m_sc_new[47]={0.0, 
                          0.0, -4.69, -20.56, -20.56, -21.65, -22.09, -38.43, -17.65, -48.34, 3.55, 
                          -38.79, -54.57, -22.76, 0.62, -54.81, -21.17, -50.51, -83.11, -45.08, -141.46, 0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
        for (int i=0; i<47; i++) sasa_m_sc[i] = sasa_m_sc_new[i];
      
        sasa_m_bb = -39.00;
        scale_factor = 10.5; // This is the found best value.
 //       scale_factor = 13.0; 
   }
   else if (solvent_type == 2){
      
        double sasa_m_sc_new[47]={0.0, 
                          0.0, -7.20, -18.75, -18.75, -41.77, -18.75, -68.10, -50.86, -85.86, -102.03, 
                          -102.03, -75.99, -67.95, -56.90, -56.90, 42.34, -65.00, -124.57, -123.19, -196.25, 
                         0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

        double sasa_b_sc_new[47]={0.0,
                        0.0, -2.28, -31.25, -31.25, -23.41, -31.25, -37.93, -27.76, -42.76, -65.73,
                        -65.73, -41.42, 0.0, -57.07, -57.07, 0.00, -85.00, -61.12, -78.71, -138.75, 
                        0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
        for (int i = 0; i <47; i++) sasa_m_sc[i] = sasa_m_sc_new[i];
        for (int i = 0; i <47; i++) sasa_b_sc[i] = sasa_b_sc_new[i];

        sasa_m_bb = -39.21;
        sasa_b_bb = -31.86;
        scale_factor =12.0;  // This is the found best value.
 //       scale_factor = 15.0;
        scale_factor_b = 1.0; 
   }
/*
{"ABC",                "GLY","ALA","SER","CYS","VAL","THR","ILE","PRO","MET","ASP",
                       "ASN","LEU","LYS","GLU","GLN","ARG","HIS","PHE","TYR","TRP",

 "CYT","GUA","ADE","THY","HOH","SOD","CLA","LIP","CHO","EAM","POT","TRE", "GOL","F00", "F01","F02","F03","F04","F05","F06","F07","F08","F09","F10","F11","MEO"}
*/
   else if (solvent_type == 3){
        double sasa_m_sc_new[47]={0.0,
                          0.0, -14.64, -39.04, -39.04, -1.02, 3.57, -25.43, -137.73, -7.65, -66.67,
                          55.69, 11.62, -110.23, -83.25, 41.41, -109.27, 42.07, -9.32, -114.32, -152.87,
                         0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

        for (int i = 0; i <47; i++) sasa_m_sc[i] = sasa_m_sc_new[i];

        sasa_m_bb = 90.0;
        scale_factor = 11.0;

   }
   else if (solvent_type == 4){
//For membrane
//
//{"ABC",                "GLY","ALA","SER","CYS","VAL","THR","ILE","PRO","MET","ASP",
//                       "ASN","LEU","LYS","GLU","GLN","ARG","HIS","PHE","TYR","TRP",
        double sasa_m_sc_new[47]={0.0,
                          54.66, 0.00, 55.00, 15.28, -23.96, 28.79, -28.82, 30.03, -15.77, -74.31,
                          67.86, -41.49, 69.67, 53.61, 51.24, 55.16, 49.69, -31.55, -0.18, -2.24,
                         0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

        for (int i = 0; i <47; i++) sasa_m_sc[i] = sasa_m_sc_new[i] - 48.28;

        sasa_m_bb = 0.0;
        scale_factor = 1.0;
   }
   else if (solvent_type == 5){
//For membrane
//
//{"ABC",                "GLY","ALA","SER","CYS","VAL","THR","ILE","PRO","MET","ASP",
//                       "ASN","LEU","LYS","GLU","GLN","ARG","HIS","PHE","TYR","TRP",
        double sasa_m_sc_new[47]={0.0,
                          -87.17, -48.28, -23.10, -47.63, -55.22, -42.75, -66.00, -13.47, -10.90, 157.85,
                          85.10, -67.98, 104.79, 123.85, 42.16, 145.48, 61.02, -43.01, 18.06, 4.66,
                         0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

        for (int i = 0; i <47; i++) sasa_m_sc[i] = sasa_m_sc_new[i];

        sasa_m_bb = 0.0;
        scale_factor = 1.0;
   }   
//For protein in air
   else if (solvent_type == 6){
        double sasa_m_sc_new[47]={0.0,
                          0.0, -54.16, 108.58, 27.86, -41.46, 106.45, -41.86, 129.82, 28.82, 212.82,
                          208.55, -44.41, 176.15, 192.66, 189.23, 169.10, 204.13, 14.27, 112.57, 101.50,
                         0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

        for (int i = 0; i <47; i++) sasa_m_sc[i] = sasa_m_sc_new[i];

        sasa_m_bb = 0.0;
        scale_factor = 1.0;
   }
   else if (solvent_type == 7){
//For membrane, Fleming
//
//{"ABC",                "GLY","ALA","SER","CYS","VAL","THR","ILE","PRO","MET","ASP",
//                       "ASN","LEU","LYS","GLU","GLN","ARG","HIS","PHE","TYR","TRP",
        double sasa_m_sc_new[47]={0.0,
                       //   7.39, -38.48, 6.19, -24.27, -48.76, 4.58, -60.75, -68.22, -45.36, 30.34,
                       //   41.15, -64.67, 71.74, 1.44, 29.05, 38.58, 63.41, -70.80, -49.01, -33.66,
                          627.9, -6572, 1088.36, -4520.88, -9795.24, 879.06, -13060.32, -12934.74, -9753.38, 5776.68,
                          7995.26, -13897.52, 15990.52, 293.02, 6027.84, 8958.04, 13353.34, -15781.22, -11134.76, -8162.7,
                         0.0,
                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

        for (int i = 0; i <47; i++) sasa_m_sc[i] = sasa_m_sc_new[i];

        sasa_m_bb = 0.0;
        scale_factor = 1.0;
   }   

//   printf("solvent_type = %d \n", solvent_type); 
/*  bgo    = (struct bgo_struct*) calloc(sim.NB,sizeof(struct bgo_struct));
                if(bgo == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for bgo\n"); exit(11);}*/
//  int boxnum;
  for(int k=0; k<sim.NB; k++){
/*          bgo[k].vdw_radii  = (double*)calloc(box[k].boxns,sizeof(double));
                if(bgo[k].vdw_radii == NULL){fprintf(stdout, "ERROR: cannot allocate memory for bgo[k].vdw_radii\n"); exit(11);}
*/
          int rescount = 0;
          int atcount = 0;
          for(int m=0; m<sim.NC; m++) {
                for(int i=0; i<bp[k][m].nbox; i++) {
                  for (int kk=0; kk<mol[m].Nres; kk++) {
                        for(int j=0; j<residue[k][rescount].Nsite; j++) {
                                  sasa[k][atcount].p = p_list[residue[k][rescount].type];
                                  sasa[k][atcount].r = vdw[residue[k][rescount].type];
                                  sasa[k][atcount].r_m = vdw_measured[residue[k][rescount].type]; 
                                  sasa[k][atcount].S = 4* PI * (sasa[k][atcount].r + RPROBE) * (sasa[k][atcount].r + RPROBE); 
                                  sasa[k][atcount].A = 0.0;    
                                  tfe[k][atcount].S_SC = sasa_ref_sc[residue[k][rescount].type]; 
                                  tfe[k][atcount].S_BB = sasa_ref_bb[residue[k][rescount].type]; 
                                  tfe[k][atcount].m_SC = sasa_m_sc[residue[k][rescount].type]/1000.0; 
                                  tfe[k][atcount].m_BB = sasa_m_bb/1000.0; 
                                  tfe[k][atcount].b_SC = sasa_b_sc[residue[k][rescount].type]/1000.0; 
                                  tfe[k][atcount].b_BB = sasa_b_bb/1000.0; 
//				  printf("m_SC is %f\n",tfe[k][atcount].m_SC);
                  		  atcount++;
              		}//for j
                	rescount++;
                  }//for kk
                }//for i
          }//for m
          #ifdef WALL
          #ifdef WWALL
          if (wall[k].hydro_index[0][0] == 4.5){
             mBBbw = 0.0; 
             mSCbw = sasa_m_sc[12]/1000.0; 
             bBBbw = 0.0; 
             bSCbw = sasa_b_sc[12]/1000.0; 
             sBBbw = sasa_ref_bb[12]; 
             sSCbw = sasa_ref_sc[12]; 
             pbw = p_list[12];
          } else if (wall[k].hydro_index[0][0] == 1.5){
             mBBbw = 0.0; 
             mSCbw = sasa_m_sc[14]/1000.0; 
             bBBbw = 0.0; 
             bSCbw = sasa_b_sc[14]/1000.0; 
             sBBbw = sasa_ref_bb[14]; 
             sSCbw = sasa_ref_sc[14]; 
             pbw = p_list[14];
          }
          #endif
          #endif
  }//for k
//  fclose(hi);
}

#endif
#endif
