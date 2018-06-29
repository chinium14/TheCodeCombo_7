#ifdef TWHAM
/* ======================================================================== */
/* twham.cpp                                                                */
/*   This file contains multiple subroutines.                               */
/* twham_init                                                               */
/*   Controls the initialization of the histograms for each of the values   */
/*   in order#.output.                                                      */
/* twham_accumulate                                                         */
/*   Driver to accumulate each of the histograms.                           */
/* twham_average                                                            */
/*   Driver to calculate the average of the accumulated properties.         */
/* twham_write                                                              */
/*   Driver to write the histograms to file.                                */
/*                                                                          */
/* Written by Thomas A. Knotts IV, 22 Sep 2005.                             */
/* ======================================================================== */

#include "defines.h"

void twham_init(void){
  double bin_width = 2.0;

  /* =============================================== */
  /* Determine the energy range over all boxes.      */
  /* =============================================== */
  double xmin=ELIMIT, xmax=-ELIMIT;
  #ifdef MPI
  MPI_Status status;
  int root=0;
  if(mpi.my_rank == root){
    double *max_temp;
    double *min_temp;
     
    max_temp = (double*)calloc(mpi.p,sizeof(double));
    min_temp = (double*)calloc(mpi.p,sizeof(double));
    max_temp[mpi.my_rank]=pe_max[0];
    min_temp[mpi.my_rank]=pe_min[0];
    for(int k=1; k<mpi.p; k++){
      MPI_Recv(&max_temp[k],1,MPI_DOUBLE,k,1,MPI_COMM_WORLD,&status);
    }
    for(int k=1; k<mpi.p; k++){
      MPI_Recv(&min_temp[k],1,MPI_DOUBLE,k,2,MPI_COMM_WORLD,&status);
    }
    for(int k=0; k<mpi.p; k++){
      if(min_temp[k] < xmin) xmin = min_temp[k];
      if(max_temp[k] > xmax) xmax = max_temp[k];
    }
    xmin = xmin-50.0;
    xmax = xmax+50.0;
    free(max_temp);
    free(min_temp);
  }
  else{
    MPI_Send(&pe_max[0],1,MPI_DOUBLE,root,1,MPI_COMM_WORLD);
    MPI_Send(&pe_min[0],1,MPI_DOUBLE,root,2,MPI_COMM_WORLD);
  }
  MPI_Bcast(&xmin,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&xmax,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  #else
  for(int k=0; k<sim.NB; k++){
    if(pe_min[k] < xmin) xmin = pe_min[k];
    if(pe_max[k] > xmax) xmax = pe_max[k];

  }
  xmin = xmin-50.0;
  xmax = xmax+50.0;
  #endif
  /* =============================================== */
  /* Allocate the historams for each box with the    */
  /* range [xmin,xmax).                              */
  /* =============================================== */
  int n = (int)((xmax-xmin)/bin_width);
  for(int k=0; k<sim.NB; k++){
    h_ebond[k] = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_ebend[k] = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_eurey[k] = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_etors[k] = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_eimpr[k] = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_elj[k]   = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_eqq[k]   = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_pe[k]    = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_ke[k]    = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_etot[k]  = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_hel[k]   = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_dnc[k]   = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_con[k]   = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_con_2[k] = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_x1[k]    = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_x2[k]    = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_rmsd[k]  = tak_histogram_calloc_uniform(n,xmin,xmax);
    h_gyr[k]   = tak_histogram_calloc_uniform(n,xmin,xmax);
  }
}
static void twham_file_write(int ibox, char file_name[], tak_histogram *h){
  const int k = ibox;
  FILE *fp;
  char path1[80];
  char path2[80];
  sprintf(path1,"./OUTPUT/BOX%d/TWHAM/",k);
  sprintf(path2,"_hist%d.dat",k);
  strcat(strcat(path1,file_name),path2);
  fp = fopen(path1, "w");
  tak_histogram_fwrite(fp,h);
  fclose(fp);
}

void twham_write(int ibox, int flag){
  #ifdef MPI
  const int k_rank = mpi.my_rank;
  #else
  const int k_rank = ibox;
  #endif
  const int k = ibox;

    if(flag == 0){
    const double max_io_time_out = 1800;
    time_t now;
    time(&now);
 if(difftime(now,wham_time[k]) < max_io_time_out) return;
 }
 
 twham_file_write(k_rank,"pe",h_pe[k]);
  //tak_histogram *h_temp = tak_histogram_div(h_ebond[k],h_pe[k]); 
  //twham_file_write(k_rank,"ebond",h_temp);
  //tak_histogram_free(h_temp);

  //tak_histogram *h_temp1 = tak_histogram_div(h_ebend[k],h_pe[k]);
  //twham_file_write(k_rank,"ebend",h_temp1);
  //tak_histogram_free(h_temp1);

  //tak_histogram *h_temp2 = tak_histogram_div(h_eurey[k],h_pe[k]);
  //twham_file_write(k_rank,"urey",h_temp2);
  //tak_histogram_free(h_temp2);

  //tak_histogram *h_temp3 = tak_histogram_div(h_etors[k],h_pe[k]);
  //twham_file_write(k_rank,"etors",h_temp3);
  //tak_histogram_free(h_temp3);

  //tak_histogram *h_temp4 = tak_histogram_div(h_eimpr[k],h_pe[k]);
  //twham_file_write(k_rank,"eimpr",h_temp4);
  //tak_histogram_free(h_temp4);
  
  //tak_histogram *h_temp5 = tak_histogram_div(h_elj[k],h_pe[k]);
  //twham_file_write(k_rank,"elj"  ,h_temp5);
  //tak_histogram_free(h_temp5);

  //tak_histogram *h_temp6 = tak_histogram_div(h_eqq[k],h_pe[k]); 
  //twham_file_write(k_rank,"eqq"  ,h_temp6);
  //tak_histogram_free(h_temp6);

  //tak_histogram *h_temp7 = tak_histogram_div(h_ke[k],h_pe[k]); 
  //twham_file_write(k_rank,"ke"   ,h_temp7);
  //tak_histogram_free(h_temp7);

  //tak_histogram *h_temp8 = tak_histogram_div(h_etot[k],h_pe[k]); 
  //twham_file_write(k_rank,"etot" ,h_temp8);
  //tak_histogram_free(h_temp8);

  //tak_histogram *h_temp9 = tak_histogram_div(h_hel[k],h_pe[k]); 
  //twham_file_write(k_rank,"hel"  ,h_temp9);
  //tak_histogram_free(h_temp9);

  tak_histogram *h_temp10 = tak_histogram_div(h_dnc[k],h_pe[k]); 
  twham_file_write(k_rank,"dnc"  ,h_temp10);
  tak_histogram_free(h_temp10);

  tak_histogram *h_temp11 = tak_histogram_div(h_con[k],h_pe[k]); 
  twham_file_write(k_rank,"con"  ,h_temp11);
  tak_histogram_free(h_temp11);

  tak_histogram *h_temp12 = tak_histogram_div(h_con_2[k],h_pe[k]); 
  twham_file_write(k_rank,"con_2",h_temp12);
  tak_histogram_free(h_temp12);

  tak_histogram *h_temp13 = tak_histogram_div(h_x1[k],h_pe[k]); 
  twham_file_write(k_rank,"x1"   ,h_temp13);
  tak_histogram_free(h_temp13);

  //tak_histogram *h_temp14 = tak_histogram_div(h_x2[k],h_pe[k]); 
  //twham_file_write(k_rank,"x2"   ,h_temp14);
  //tak_histogram_free(h_temp14);

  //tak_histogram *h_temp15 = tak_histogram_div(h_rmsd[k],h_pe[k]); 
  //twham_file_write(k_rank,"rmsd" ,h_temp15);
  //tak_histogram_free(h_temp15);

  tak_histogram *h_temp16 = tak_histogram_div(h_gyr[k],h_pe[k]); 
  twham_file_write(k_rank,"gry"  ,h_temp16);
  tak_histogram_free(h_temp16);
  /*twham_file_write(k_rank,"ebond",tak_histogram_div(h_ebond[k],h_pe[k]));
  twham_file_write(k_rank,"ebend",tak_histogram_div(h_ebend[k],h_pe[k]));
  twham_file_write(k_rank,"eurey",tak_histogram_div(h_eurey[k],h_pe[k]));
  twham_file_write(k_rank,"etors",tak_histogram_div(h_etors[k],h_pe[k]));
  twham_file_write(k_rank,"eimpr",tak_histogram_div(h_eimpr[k],h_pe[k]));
  twham_file_write(k_rank,"elj"  ,tak_histogram_div(h_elj[k]  ,h_pe[k]));
  twham_file_write(k_rank,"eqq"  ,tak_histogram_div(h_eqq[k]  ,h_pe[k]));
  twham_file_write(k_rank,"ke"   ,tak_histogram_div(h_ke[k]   ,h_pe[k]));
  twham_file_write(k_rank,"etot" ,tak_histogram_div(h_etot[k] ,h_pe[k]));
  twham_file_write(k_rank,"hel"  ,tak_histogram_div(h_hel[k]  ,h_pe[k]));
  twham_file_write(k_rank,"dnc"  ,tak_histogram_div(h_dnc[k]  ,h_pe[k]));
  twham_file_write(k_rank,"con"  ,tak_histogram_div(h_con[k]  ,h_pe[k]));
  twham_file_write(k_rank,"con_2",tak_histogram_div(h_con_2[k],h_pe[k]));
  twham_file_write(k_rank,"x1"   ,tak_histogram_div(h_x1[k]   ,h_pe[k]));
  twham_file_write(k_rank,"x2"   ,tak_histogram_div(h_x2[k]   ,h_pe[k]));
  twham_file_write(k_rank,"rmsd" ,tak_histogram_div(h_rmsd[k] ,h_pe[k]));
  twham_file_write(k_rank,"gry"  ,tak_histogram_div(h_gyr[k]  ,h_pe[k])); */
  time(&wham_time[k]);
}
void twham_accumulate(int ibox){
  const int k = ibox;
  const double pe = en[k].potens;
  struct energy *enk = &en[k];
  struct struct_quant *ordk = &ordparam[k];
  tak_histogram_accumulate(h_ebond[k],pe,enk->bond);
  tak_histogram_accumulate(h_ebend[k],pe,enk->bend);
  #ifndef NEUTRAL
  tak_histogram_accumulate(h_eurey[k],pe,enk->urey);
  #endif
  tak_histogram_accumulate(h_etors[k],pe,enk->tors);
  tak_histogram_accumulate(h_eimpr[k],pe,enk->impr);
  tak_histogram_accumulate(h_elj[k]  ,pe,(enk->nbonds-enk->coulomb));
  tak_histogram_accumulate(h_eqq[k]  ,pe,enk->coulomb);
  #ifndef SASAREX
  tak_histogram_increment(h_pe[k]   ,pe);
  #else
  tak_histogram_increment(h_pe[k]   ,(enk->esasa));
  #endif
  tak_histogram_accumulate(h_ke[k]   ,pe,enk->kinet);
  tak_histogram_accumulate(h_etot[k] ,pe,enk->totals);
  tak_histogram_accumulate(h_hel[k]  ,pe,ordk->hel);
  tak_histogram_accumulate(h_dnc[k]  ,pe,ordk->d_nc);  
  #ifndef SASAREX
  tak_histogram_accumulate(h_con[k]  ,pe,ordk->con);
  #else
  tak_histogram_accumulate(h_con[k]  ,enk->esasa, ordk->con);
  #endif
  tak_histogram_accumulate(h_con_2[k],pe,ordk->con_2); 
  tak_histogram_accumulate(h_x1[k]   ,pe,ordk->x1);
  tak_histogram_accumulate(h_x2[k]   ,pe,ordk->x2);
  tak_histogram_accumulate(h_rmsd[k] ,pe,ordk->rmsd);
  tak_histogram_accumulate(h_gyr[k]  ,pe,ordk->gyr);



}



#endif
