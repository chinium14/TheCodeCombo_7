#ifdef CUDA
#ifndef PRESSURE
#include "cuda.h"

__global__ void nonbond_loop(double * force, double * enonbond, double * enonbonds, int boxns, real3 * atom_xyz, double bl, double bh, interb *ljset, int k_dev){
	int a = threadIdx.x + blockIdx.x * blockDim.x;
	int b = threadIdx.y + blockIdx.y * blockDim.y;
	int id = a + b * blockDim.x * gridDim.x;
	if (b >= a){
		double fxa = 0.0;
		double fya = 0.0;
		double fza = 0.0;
		double ax = atom_xyz[a].x;
		double ay = atom_xyz[a].y;
		double az = atom_xyz[a].z;
		double bx = atom_xyz[b].x;
		double by = atom_xyz[b].y;
		double bz = atom_xyz[b].z;
		double dx = ax - bx;
		double dy = ay - by;
		double dz = az - bz;
		double pbcx = 0.0;
		double pbcy = 0.0;
		double pbcz = 0.0;
                
		if(dx >  bh) pbcx =- bl;
                if(dx < -bh) pbcx =+ bl;
                if(dy >  bh) pbcy =- bl;
                if(dy < -bh) pbcy =+ bl;
                if(dz >  bh) pbcz =- bl;
                if(dz < -bh) pbcz =+ bl;
		
		dx += pbcx;
                dy += pbcy;
                dz += pbcz;
                double dr2 = dx*dx + dy*dy + dz*dz;
                if(dr2 < rc2) {
        
//      //      printf("%d %d  %lf \n",a,b,sqrt(dr2));
		double eps = ljset[k].pot[a][b].eps;
                double sig = ljset[k].pot[a][b].sig;
                double qq  = ljset[k].pot[a][b].qq;
        
                double sr2  = (sig * sig) / dr2;
                double sr6  = sr2 * sr2 * sr2;
                double sr12 = sr6 * sr6;
        
                double sr2s  = (sig * sig) / rc2;
                double sr6s  = sr2s * sr2s * sr2s;
                double sr12s = sr6s * sr6s;
		
		force[id] = (12.0*eps/dr2) * (sr12 -sr6);
		enonbond[id] += eps * (sr12 - 2*sr6);
		enonbonds[id] += eps * (sr12 - 2*sr6);
		}else{
			force[id] = 0;
			enonbond[id] = 0;
			enonnbonds[id] = 0;
		}
		__syncthreads();
		int half = boxns/2;
		while (half != 0){
			force[id] +=force[id+half];
			__syncthreads();
			half /= 2;
		}
	}
}

void cnonbond(int ibox, double *energy)
{
	int k = ibox;
	double enonbond = 0.0;
	double enonbonds = 0.0;
	double ecoulomb = 0.0;
	
	double enonbond_temp;
	double enonbonds_temp;
	double ecoulomb_temp;
	
	int k_dev, boxns_dev;
	double *enonbond_dev;
	double *enonbonds_dev;
	double *force_dev;
	cudaMalloc((void*)&k_dev,sizeof(int));
	cudaMalloc((void*)&boxns_dev,sizeof(int));
	cudaMalloc((void**)&enonbond_dev,boxns[k].boxns*boxns[k].boxns*sizeof(double));
	cudaMalloc((void**)&enonbonds_dev,boxns[k].boxns*boxns[k].boxns*sizeof(double));
	cudaMalloc((void**)&force_dev,boxns[k].boxns*boxns[k].boxns*sizeof(double));
	cudaMemcpy(k_dev,k,sizeof(int),cudaMemcpyHostToDevice);	
	cudaMemcpy(boxns_dev,box[k].boxns,sizeof(int),cudaMemcpyHostToDevice);
	
	  double rc   = sqrt(box[k].rc2);
	  double rc2  = box[k].rc2;
	  double rfcs = (sim.epsRF[k] - 1.0) / (2.0*sim.epsRF[k] + 1.0);
	  double rfc  = rfcs / (rc2*rc);
	  double bl   = box[k].boxl;
	  double bh   = box[k].boxh;
		
?	  double bl_dev;
?	  double bh_dev;
	  cudaMemcpy(bl_dev,bl, sizeof(double),cudaMemcpyHostToDevice);
	  cudaMemcpy(bh_dev,bh, sizeof(double),cudaMemcpyHostToDevice);
//Generate atom coordinate information
	real3 * atom_xyz;
	atom_xyz = (real3 *)malloc(box[k].boxns * sizeof(real3));
	for (int i = 0; i < box[k].boxns-1; i++){
		atom_xyz[i] = make_real3(atom[k][i].x, atom[k][i].y, atom[k][i].z)
	}
	real3 * atom_xyz_dev;
	cudaMalloc((void **)&atom_xyz_dev, box[k].boxns*sizeof(real3));
	cudaMemcpy(atom_xyz_dev,atom_xyz,box[k].boxns*sizeof(real3),cudaMemcpyHostToDevice);
// Perform calculation on GPU of nonbond forces
	dim3 threads(16,16);
	dim3 blocks((boxns+15)/16,(boxns+15)/16);
	nonbond_loop<<<blocks,threads>>>(force_dev,enonbond_dev,enonbonds_dev,boxns,atom_xyz_dev, bl_dev, bh_dev, );

	cudaMemcpy(k,k_dev,sizeof(int),cudaMemcpyDeviceToHost);	
	cudaMemcpy(enonbond,enonbond_dev[0],sizeof(double),cudaMemcpyDeviceToHost);	
	cudaMemcpy(enonbonds,enonbonds_dev[0],sizeof(double),cudaMemcpyDeviceToHost);	
	cudaMemcpy(force,force_dev[0],sizeof(double),cudaMemcpyDeviceToHost);
	cudaFree(k_dev);
	cudaFree(enonbond_dev);
	cudaFree(enonbonds_dev);
	cudaFree(force_dev);
	cudaFree(atom_xyz_dev);
	free(atom_xyz);
}
#endif
#endif
