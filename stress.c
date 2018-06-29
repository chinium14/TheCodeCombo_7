/********************************************************/
/* stress.c                                             */
/*                                                      */
/* Contains routines for calculating stress invariants, */
/* computing block averages and errors, and outputting  */
/* data. Contains 6 functions.                          */
/*                                                      */
/*======================================================*/
/*    stressComputeInvariantsNow(int ibox)              */
/*                                                      */
/*       Input:      Box Number                         */
/*                                                      */
/*       Output:     Stress invariants (I1, I2, I3)     */
/*                   Mean stress (P)                    */
/*                   Eigenvalues of stress tensor       */
/*                   Deviatoric invariants (J1, J2, J3) */
/*                   Eigenvalues of deviatoric tensor   */
/*                   von Mises stress (VMS)             */
/*                                                      */
/*======================================================*/
/*    stressSumInvariants(int ibox)                     */
/*                                                      */
/*       Input:      Box number                         */
/*                                                      */
/*       Accumulates I_sum. Called each time step.      */
/*                                                      */
/*======================================================*/
/*    stressSumTensors(int ibox)                        */
/*                                                      */
/*       Input:      Box number                         */
/*                                                      */
/*       Accumulates S_sum. Called each time step.      */
/*                                                      */
/*======================================================*/
/*    stressOutputInvariantsNow(int ibox)               */
/*                                                      */
/*       Input:      Box number                         */
/*                                                      */
/*       Outputs instantaneous stress invariants for a  */
/*       single time step. Mainly for error checking.   */
/*                                                      */
/*======================================================*/
/*    stressOutputTensorsNow(int ibox)                  */
/*                                                      */
/*       Input:      Box number                         */
/*                                                      */
/*       Outputs instantaneous stress tensors for a     */
/*       single time step. Mainly for error checking.   */
/*                                                      */
/*======================================================*/
/*    stressFinishBlock(int ibox)                       */
/*                                                      */
/*       Input:      Box number                         */
/*                                                      */
/*       Called at the end of each block. Keeps running */
/*       totals of I_bar_sum and I_bar_sumsq.           */
/*                                                      */
/*======================================================*/
/*    stressResetBlock(int ibox)                        */
/*                                                      */
/*       Input:      Box number                         */
/*                                                      */
/*       Called after stressFinishBlock(). Resets I_sum */
/*       to zero.                                       */
/*                                                      */
/*======================================================*/
/*    stressOutputStats(int ibox)                       */
/*                                                      */
/*       Input:      Box number                         */
/*                                                      */
/*       Outputs averages and errors at the end of the  */
/*       production run.                                */
/*                                                      */
/********************************************************/
     
#ifdef ZHOU

#include "defines.h"
//#include "eig3.h"

double rms (int, double, double);

void stressComputeInvariantsNow(int ibox)
{
//	printf("Computing invariants now\n");
        int k = ibox;

        for(int aa = 0; aa<box[k].boxns;aa++)
        {
                double I1,I2,I3,J1,J2,J3,P;
                double Eval_I[3];

		stresses[k][aa][0][1] = stresses[k][aa][1][0];
		stresses[k][aa][0][2] = stresses[k][aa][2][0];
		stresses[k][aa][1][2] = stresses[k][aa][2][1];

                double s11 = stresses[k][aa][0][0], s12 = stresses[k][aa][0][1], s13 = stresses[k][aa][0][2];
                double s21 = stresses[k][aa][1][0], s22 = stresses[k][aa][1][1], s23 = stresses[k][aa][1][2];
                double s31 = stresses[k][aa][2][0], s32 = stresses[k][aa][2][1], s33 = stresses[k][aa][2][2];

                //Do the things to compute stress invariants.
                I1 = s11+s22+s33;       //Trace
                I2 = s11*s22 + s22*s33 + s33*s11 - s12*s21 - s23*s32 - s31*s13;
                I3 = s11*(s22*s33-s23*s32) - s12*(s21*s33-s23*s31) + s13*(s21*s32-s22*s31);  //Determinant
                P =  I1/3.0;  //Mean stress = 1/3 Trace
//                J1 = 0.0;  //Trace = 0
//                J2 = (s11-P)*(s22-P) + (s22-P)*(s33-P) + (s33-P)*(s11-P) - s12*s21 - s23*s32 - s31*s13;
//                J3 = (s11-P)*((s22-P)*(s33-P)-s23*s32) - s12*(s21*(s33-P)-s23*s31) + s13*(s21*s32-(s22-P)*s31);  //Determinant

		//Shortcuts
		J1 = 0.0;
		J2 = I1*I1/3.0 - I2;
		J3 = I1*I1*I1/13.5 - I1*I2/3.0 + I3;

                //Do eigenvalues:

                //Characteristic polynomial for stress tensor:
                //      -E^3 + I1 E^2 - I2 E + I3 = 0
                //       x^3 +  a x^2 +  b x +  c = 0

                //Following procedure for finding the roots of a cubic taken from Numerical Recipes:

                int cpxI=0,cpxJ=0;
                double a,b,c,Q,R,Q3,q2,x1=0,x2=0,x3=0;
                a = -I1;
                b = I2;
                c = -I3;
                Q = (a*a-3.0*b)/9.0;
                R = (2.0*a*a*a-9.0*a*b+27.0*c)/54.0;
                Q3 = Q*Q*Q;
                q2 = sqrt(Q);
                x1=0,x2=0,x3=0;
                if(R*R<Q3)
                {
                        //Three real solutions
                        cpxI = 0;
                        double theta=acos(R/sqrt(Q3));
                        x1 = -2.0*q2*cos( theta      /3.0)-a/3.0;
                        x2 = -2.0*q2*cos((theta+2*PI)/3.0)-a/3.0;
                        x3 = -2.0*q2*cos((theta-2*PI)/3.0)-a/3.0;
                        Eval_I[0]=x1;
                        Eval_I[1]=x2;
                        Eval_I[2]=x3;
                }else{
                      //One real, two complex solutions
                        cpxI = 1;
                        double A;
                        if(R>0)
                                A = -pow(fabs(R)+sqrt(R*R-Q3),0.3333333333333333);
                        else
                                A = pow(fabs(R)+sqrt(R*R-Q3),0.3333333333333333);
                                double B;
                        if(A==0) B=0; else B=Q/A;
                        x1 = (A+B)-a/3.0;               //Real solution
                        x2 = -0.5*(A+B)-a/3.0;          //Real part of cmplx solutions
                        x3 = sqrt(3.0)/2.0*(A-B);       //Imag part of cmplx solutions
                        Eval_I[0]=x1;
                        Eval_I[1]=x2;
                        Eval_I[2]=x3;
                }
/*
                x1=0,x2=0,x3=0;
                a = -0.0;
                b = J2;
                c = -J3;
                Q = (a*a-3.0*b)/9.0;
                R = (2.0*a*a*a-9.0*a*b+27.0*c)/54.0;
                Q3 = Q*Q*Q;
                q2 = sqrt(Q);
                x1=0,x2=0,x3=0;
                if(R*R<Q3)
                {
                        //Three real solutions
                        cpxJ = 0;
                        double theta=acos(R/sqrt(Q3));
                        x1 = -2.0*q2*cos( theta      /3.0)-a/3.0;
                        x2 = -2.0*q2*cos((theta+2*PI)/3.0)-a/3.0;
                        x3 = -2.0*q2*cos((theta-2*PI)/3.0)-a/3.0;
                        Eval_J[0]=x1;
                        Eval_J[1]=x2;
                        Eval_J[2]=x3;
                }else{
                        //One real, two complex solutions
                        cpxJ = 1;
                        double A;
                        if(R>0)
                                A = -pow(fabs(R)+sqrt(R*R-Q3),0.3333333333333333);
                        else
                                A = pow(fabs(R)+sqrt(R*R-Q3),0.3333333333333333);
                        double B;
                        if(A==0) B=0; else B=Q/A;
                        x1 = (A+B)-a/3.0;         //Real solution
                        x2 = -0.5*(A+B)-a/3.0;    //Real part of complex solutions
                        x3 = sqrt(3.0)/2.0*(A-B); //Imag part of complex solutions
                        Eval_J[0]=x1;
                        Eval_J[1]=x2;
                        Eval_J[2]=x3;
                }

*/

                stress_i[k][aa].I1   = I1;
                stress_i[k][aa].I2   = I2;
                stress_i[k][aa].I3   = I3;
                stress_i[k][aa].J1   = J1;
                stress_i[k][aa].J2   = J2;
                stress_i[k][aa].J3   = J3;
                stress_i[k][aa].P    = P;
                stress_i[k][aa].cpxI = cpxI;
                stress_i[k][aa].cpxJ = cpxJ;
                stress_i[k][aa].VMS  = sqrt(3.0 * J2);
	
		//If x is an eigenvalue of I, then x - P is an eigenvalue of J.		
		if(cpxI == 0)
		{	
			//Real solutions
               		for(int i=0;i<3;i++)
	                {
	                        stress_i[k][aa].Eval_I[i] = Eval_I[i];
	                        stress_i[k][aa].Eval_J[i] = Eval_I[i] - P;	
			}
		}
		else
		{
			//Two complex, one real solution
			for(int i=0;i<2;i++)
			{
	                        stress_i[k][aa].Eval_I[i] = Eval_I[i];		
	                        stress_i[k][aa].Eval_J[i] = Eval_I[i] - P;	//This is the real part of the complex solutions, or the real solution
			}			
	                stress_i[k][aa].Eval_I[2] = Eval_I[2];
	                stress_i[k][aa].Eval_J[2] = Eval_I[2];		//No "- P" is intentional because this is the imaginary part of the complex solutions
		}

        } //Loop over all sites.

        return;
}

void stressSumInvariants(int ibox)
{
//	printf("Summing invariants\n");
	int k = ibox;
	for(int aa=0;aa<box[k].boxns;aa++)
	{
		stress_isum[k][aa].I1  += stress_i[k][aa].I1;
		stress_isum[k][aa].I2  += stress_i[k][aa].I2;
		stress_isum[k][aa].I3  += stress_i[k][aa].I3;
		stress_isum[k][aa].P   += stress_i[k][aa].P;
		stress_isum[k][aa].J1  += stress_i[k][aa].J1;
		stress_isum[k][aa].J2  += stress_i[k][aa].J2;
		stress_isum[k][aa].J3  += stress_i[k][aa].J3;
	 	//Not sure how to sum eigenvalues yet...
	} //Loop over all sites in box k
	return;
}

void stressSumTensors(int ibox)
{
//	printf("Summing tensors\n");
	int k = ibox;
	for(int aa = 0; aa<box[k].boxns; aa++)
		for(int x=0;x<3;x++)
			for(int y=0;y<3;y++)
			{
				stresses[k][aa][0][1] = stresses[k][aa][1][0];
				stresses[k][aa][0][2] = stresses[k][aa][2][0];
				stresses[k][aa][1][2] = stresses[k][aa][2][1];
				stresses_sum[k][aa][x][y] += stresses[k][aa][x][y];
			}
	//Loop over all sites in box k

	return;
}

void stressOutputInvariantsNow(int ibox, unsigned int icycle, int flag)
{
//	printf("Outputting invariants\n");
	int k = ibox;
	
	char namefile[80];
        FILE *fid;
	if (stressOutputCount< 10)
                sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stress000%d.inv",k,stressOutputCount);
        else if (stressOutputCount < 100)
                sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stress00%d.inv",k,stressOutputCount);
        else if (stressOutputCount < 1000)
                sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stress0%d.inv",k,stressOutputCount);
        else
                sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stress%d.inv",k,stressOutputCount);

        fid = fopen(namefile,"w");


	//Calculate the simulation time.
	double tt;
	if (flag==0) tt = sim.dt * 1.0*icycle;
	else if (flag==1 && (sim.ID2==3 || sim.ID2==5)) tt= sim.dt*sim.cyc_eq+ sim.dtlong*1.0*icycle; // for multiple time step
	else tt = sim.dt * 1.0*(sim.cyc_eq+icycle);

        fprintf(fid,"Stress invariants for each site, box %d, time = %lf\n\n",k,tt);
        fprintf(fid,"Site   \tP           \tI1          \tI2          \tI3          \tJ1          \tJ2          \tJ3          \tVMS\n");

        for (int i=0;i<box[k].boxns;i++)
        {

        //            %12.7lf\t         Width = 12, precision = 7, long float, with a tab.
        //            %7lu\t            Width = 7, long int, with a tab.
                fprintf(fid,"%7i\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t",
                        i,
                        stress_i[k][i].P,
                        stress_i[k][i].I1,
                        stress_i[k][i].I2,
                        stress_i[k][i].I3,
                        stress_i[k][i].J1,
                        stress_i[k][i].J2,
                        stress_i[k][i].J3,
                        stress_i[k][i].VMS);
                fprintf(fid,"\n");

        }

        fclose(fid);

	return;
}

void stressOutputEigenvaluesNow(int ibox, unsigned int icycle, int flag)
{
//	printf("Outputting eigenvalues\n");
        int k = ibox;

        char namefile[80];
        FILE *fid;
        if (stressOutputCount< 10)
                sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stress000%d.eig",k,stressOutputCount);
        else if (stressOutputCount < 100)
                sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stress00%d.eig",k,stressOutputCount);
        else if (stressOutputCount < 1000)
                sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stress0%d.eig",k,stressOutputCount);
        else
                sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stress%d.eig",k,stressOutputCount);

        fid = fopen(namefile,"w");
	
	//Calculate the current simulation time.
	double tt;
	if (flag==0) tt = sim.dt * 1.0*icycle;
	else if (flag==1 && (sim.ID2==3 || sim.ID2==5)) tt= sim.dt*sim.cyc_eq+ sim.dtlong*1.0*icycle; // for multiple time step
	else tt = sim.dt * 1.0*(sim.cyc_eq+icycle);

        fprintf(fid,"Eigenvalues for each site, box %d, time = %lf\n",k,tt);
	fprintf(fid,"Order is not guaranteed to be consistant!\n\n");
        
	fprintf(fid,"Site \tEI1(real) \tEI1(imag) \tEI2(real) \tEI2(imag) \tEI3(real) \tEI3(imag) \t");	//Output continues....
	       fprintf(fid,"EJ1(real) \tEJ1(imag) \tEJ2(real) \tEJ2(imag) \tEJ3(real) \tEJ3(imag)\n");

        for (int i=0;i<box[k].boxns;i++)
        {

        //            %12.7lf\t         Width = 12, precision = 7, long float, with a tab.
        //            %7lu\t            Width = 7, long int, with a tab.
		if(stress_i[k][i].cpxI == 0)
		{
	                fprintf(fid,"%7i\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e",
	                        i,
	                        stress_i[k][i].Eval_I[0], 0.0,
	                        stress_i[k][i].Eval_I[1], 0.0,
	                        stress_i[k][i].Eval_I[2], 0.0,
	                        stress_i[k][i].Eval_J[0], 0.0,
	                        stress_i[k][i].Eval_J[1], 0.0,
	                        stress_i[k][i].Eval_J[2], 0.0);
		} else {
	                fprintf(fid,"%7i\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e\t%10.7e",
	                        i,
	                        stress_i[k][i].Eval_I[0], 0.0,
	                        stress_i[k][i].Eval_I[1], stress_i[k][i].Eval_I[2],
	                        stress_i[k][i].Eval_I[1], -stress_i[k][i].Eval_I[2],
	                        stress_i[k][i].Eval_J[0], 0.0,
	                        stress_i[k][i].Eval_J[1], stress_i[k][i].Eval_J[2],
	                        stress_i[k][i].Eval_J[1], -stress_i[k][i].Eval_J[2]);
			
		}
                fprintf(fid,"\n\n");

        }

        fclose(fid);

        return;
}

void stressOutputTensorsNow(int ibox, unsigned int icycle, int flag)
{
//	printf("Outputting tensors\n");
        int k = ibox;
        FILE *fid;
        char namefile[80];
        if (stressOutputCount< 10)
                sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stress000%d.tensor",k,stressOutputCount);
        else if (stressOutputCount < 100)
                sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stress00%d.tensor",k,stressOutputCount);
        else if (stressOutputCount < 1000)
                sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stress0%d.tensor",k,stressOutputCount);
        else
                sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stress%d.tensor",k,stressOutputCount);

        fid = fopen(namefile,"w");

	//Calculate the current simulation time.
	double tt;
	if (flag==0) tt = sim.dt * 1.0*icycle;
	else if (flag==1 && (sim.ID2==3 || sim.ID2==5)) tt= sim.dt*sim.cyc_eq+ sim.dtlong*1.0*icycle; // for multiple time step
	else tt = sim.dt * 1.0*(sim.cyc_eq+icycle);

        fprintf(fid,"Stress tensors for each site, box %d, time = %lf\n\n",k,tt);


        for (int i=0;i<box[k].boxns;i++)
        {
		stresses[k][i][0][1] = stresses[k][i][1][0];
		stresses[k][i][0][2] = stresses[k][i][2][0];
		stresses[k][i][1][2] = stresses[k][i][2][1];
                fprintf(fid,"Site %u\n",i);
                fprintf(fid,"%12.7e\t%12.7e\t%12.7e\n",stresses[k][i][0][0],stresses[k][i][1][0],stresses[k][i][2][0]);
                fprintf(fid,"%12.7e\t%12.7e\t%12.7e\n",stresses[k][i][0][1],stresses[k][i][1][1],stresses[k][i][2][1]);
                fprintf(fid,"%12.7e\t%12.7e\t%12.7e\n",stresses[k][i][0][2],stresses[k][i][1][2],stresses[k][i][2][2]);
        }

        fclose(fid);
	return;
}

void stressFinishBlock(int ibox)
{
//	printf("Finishing block\n");
	//Gets ibarsum and ibarsumsq
	
	int k = ibox;
	int n = stressStepCount;	//  # steps in one block

	for(int aa=0;aa<box[k].boxns;aa++)
	{

		double I1bar = stress_isum[k][aa].I1 / n;
		double I2bar = stress_isum[k][aa].I2 / n;
		double I3bar = stress_isum[k][aa].I3 / n;
		double _Pbar = stress_isum[k][aa].P  / n;
		double J1bar = stress_isum[k][aa].J1 / n;
		double J2bar = stress_isum[k][aa].J2 / n;
		double J3bar = stress_isum[k][aa].J3 / n;

		stress_ibarsum[k][aa].I1 += I1bar;
		stress_ibarsum[k][aa].I2 += I2bar; 
		stress_ibarsum[k][aa].I3 += I3bar;
		stress_ibarsum[k][aa].P  += _Pbar;
		stress_ibarsum[k][aa].J1 += J1bar;
		stress_ibarsum[k][aa].J2 += J2bar;
		stress_ibarsum[k][aa].J3 += J3bar;

		stress_ibarsumsq[k][aa].I1 += I1bar * I1bar;
		stress_ibarsumsq[k][aa].I2 += I2bar * I2bar; 
		stress_ibarsumsq[k][aa].I3 += I3bar * I3bar;
		stress_ibarsumsq[k][aa].P  += _Pbar * _Pbar;
		stress_ibarsumsq[k][aa].J1 += J1bar * J1bar;
		stress_ibarsumsq[k][aa].J2 += J2bar * J2bar;
		stress_ibarsumsq[k][aa].J3 += J3bar * J3bar;
	}

	stressBlockCount++;
	return;
}

void stressResetBlock(int ibox)
{	
//	printf("Resetting block\n");
	stressStepCount = 0;

	int k = ibox;
	for(int aa=0; aa<box[k].boxns; aa++)
	{
		stress_isum[k][aa].I1 = 0;
		stress_isum[k][aa].I2 = 0;
		stress_isum[k][aa].I3 = 0;
		stress_isum[k][aa].P  = 0;
		stress_isum[k][aa].J1 = 0;
		stress_isum[k][aa].J2 = 0;
		stress_isum[k][aa].J3 = 0;
	}

	return;
}

void stressResetTensors(int ibox)
{
//	printf("Resetting tensors\n");
	int k = ibox;
	for(int aa = 0; aa<box[k].boxns; aa++)
		for(int x=0;x<3;x++)
			for(int y=0;y<3;y++)
				stresses[k][aa][x][y] = 0;
	//Loop over all sites in box k

	return;
}

void stressResetInvariant( struct stress_type inv )
{
//	printf("Resetting invariants\n");
	inv.I1 = 0;
	inv.I2 = 0;
	inv.I3 = 0;
	inv.P  = 0;
	inv.J1 = 0;
	inv.J2 = 0;
	inv.J3 = 0;

	return;
}



		
void stressResetAll(int ibox)
{
//	printf("Resetting all\n");
	//Reset everything
		
	int k = ibox;

/*	#ifdef ZEROAM 
		char namefile[80],namefile2[80];
	        FILE *fid,*fid2;
	        sprintf(namefile, "./OUTPUT/BOX%d/STRESSES/angmom.dat",k);
	        fid = fopen(namefile,"w");      //Erase the old angmom.dat file
		fclose(fid);
	        sprintf(namefile2, "./OUTPUT/BOX%d/STRESSES/linmom.dat",k);
	        fid2 = fopen(namefile2,"w");      //Erase the old linmom.dat file
		fclose(fid2);
	#endif
*/
	stressResetTensors(ibox);
	stressResetBlock(ibox);
	stressBlockCount = 0;
	stressTotalCount = 0;
	stressOutputCount = 0;
	for(int aa = 0; aa<box[k].boxns; aa++)
	{
		stressResetInvariant(stress_ibarsum  [k][aa]);
		stressResetInvariant(stress_ibarsumsq[k][aa]);

		for(int x=0;x<3;x++)  for(int y=0;y<3;y++)  stresses_sum[k][aa][x][y] = 0;

	} //Loop over all sites in box k

	return;	
}

void stressOutputTensorAverages(int ibox)
{
//	printf("Outputting tensor averages\n");
	int k = ibox;
	FILE *fid;
        char namefile[80];
        sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stressmean.tensor",k);
        fid = fopen(namefile,"w");
        fprintf(fid,"Tensor averages for each site, box %d\n",k);
	#ifndef ZEROAM
		fprintf(fid,"WARNING: Angular momentum might not have zeroed. Data will be confounded if molecule has rotated.\n");
	#endif
	fprintf(fid,"\n");
	int n = stressTotalCount;
	for (int i=0;i<box[k].boxns;i++)
        {
                fprintf(fid,"Site %u\n",i);
                fprintf(fid,"%12.7e\t%12.7e\t%12.7e\n",stresses_sum[k][i][0][0]/n,stresses_sum[k][i][1][0]/n,stresses_sum[k][i][2][0]/n);
                fprintf(fid,"%12.7e\t%12.7e\t%12.7e\n",stresses_sum[k][i][0][1]/n,stresses_sum[k][i][1][1]/n,stresses_sum[k][i][2][1]/n);
                fprintf(fid,"%12.7e\t%12.7e\t%12.7e\n",stresses_sum[k][i][0][2]/n,stresses_sum[k][i][1][2]/n,stresses_sum[k][i][2][2]/n);
        }

        fclose(fid);
        return;
}
void stressOutputStats(int ibox)
{
//	printf("Outputting stats\n");
	int k = ibox;
	
	FILE *fid;
        char namefile[80];
        sprintf(namefile,"./OUTPUT/BOX%d/STRESSES/stress.stats",k);

        fid = fopen(namefile,"w");

        fprintf(fid,"Stress invariants for each site, box %d\n\n",k);

	fprintf(fid,"Site \tI1          \terr         \tI2          \terr         \tI3          \terr         \tP           \terr         \tJ1          \terr         \tJ2          \terr         \tJ3          \terr\n");

	int b = stressBlockCount;
	for(int aa = 0; aa<box[k].boxns; aa++)
	{

		double I1 = stress_ibarsum[k][aa].I1 / b;
		double I2 = stress_ibarsum[k][aa].I2 / b;
		double I3 = stress_ibarsum[k][aa].I3 / b;
		double P  = stress_ibarsum[k][aa].P  / b;
		double J1 = stress_ibarsum[k][aa].J1 / b;
		double J2 = stress_ibarsum[k][aa].J2 / b;
		double J3 = stress_ibarsum[k][aa].J3 / b;
		double I1err = rms(b, stress_ibarsum[k][aa].I1, stress_ibarsumsq[k][aa].I1); 
		double I2err = rms(b, stress_ibarsum[k][aa].I2, stress_ibarsumsq[k][aa].I2); 
		double I3err = rms(b, stress_ibarsum[k][aa].I3, stress_ibarsumsq[k][aa].I3); 
		double _Perr = rms(b, stress_ibarsum[k][aa].P , stress_ibarsumsq[k][aa].P ); 
		double J1err = rms(b, stress_ibarsum[k][aa].J1, stress_ibarsumsq[k][aa].J1); 
		double J2err = rms(b, stress_ibarsum[k][aa].J2, stress_ibarsumsq[k][aa].J2); 
		double J3err = rms(b, stress_ibarsum[k][aa].J3, stress_ibarsumsq[k][aa].J3); 
		fprintf(fid,"%7i\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t",
			aa,
			I1,I1err,
			I2,I2err,
			I3,I3err,
			P ,_Perr,
			J1,J1err,
			J2,J2err,
			J3,J3err);
		fprintf(fid,"\n");
        }

        fclose(fid);
	stressOutputTensorAverages(k);
	return;
}

#ifdef TRR
void stressSumMovieStress()
{
	//You'll note that stressSumInvariants needs to be called before this one, for obvious reasons.
	for(int aa=0; aa<box[0].boxns; aa++)
		stress_movie_psum[aa] += stress_i[0][aa].P;
}
void stressInitMovie(int ibox)
{
	int k = ibox;
	stressFrameCount = 0;
	//Get the movie file ready for writing
	char namefile3[80];
	FILE *fid3;
	sprintf(namefile3, "./OUTPUT/BOX%d/STRESSES/stress.movie",k);
	if ( !(fid3 = fopen(namefile3,"w")) ) ;
		fclose(fid3);	//Erase the old stress.movie file

	return;	
}

void stressOutputFrame(int ibox, unsigned int icycle, int flag) 
{
	//This routine outputs instantaneous stress values for use in 
	//making movies with colors that change. Right now, it does
	//only the mean stress, and the instantaneous value.

	//As of NOW, this function will output the average accumulated mean stress, since the last frame.

	int k = ibox;

	//Open file for writing
	char namefile[80];
	FILE *fid;
	sprintf(namefile, "./OUTPUT/BOX%d/STRESSES/stress.movie",k);
	fid = fopen(namefile,"a");

	//Calculate the simulation time.
	double tt;
	if (flag==0) tt = sim.dt * 1.0*icycle;
	else if (flag==1 && (sim.ID2==3 || sim.ID2==5)) tt= sim.dt*sim.cyc_eq+ sim.dtlong*1.0*icycle; // for multiple time step
	else tt = sim.dt * 1.0*(sim.cyc_eq+icycle);

	//Output the next frame of information
        fprintf(fid,"Frame %i Time %f\n",stressFrameCount, tt);
	double iblockt = 1.0 / sim.blockt;
	for (int i=0;i<box[k].boxns;i++)
	{
		
                fprintf( fid,"%12.7e\n",stress_movie_psum[i] * iblockt );	//Output mean stress
		stress_movie_psum[i] = 0;					//Reset for next time
	}

        fclose(fid);

	stressFrameCount += 1;
	return;
}
#endif //TRR

#ifdef ZEROAM

void outputAngularMomentum(int ibox, unsigned int icycle, int flag)
{
	int k = ibox;
        char namefile[80],namefile2[80];
        FILE *fid,*fid2;
        sprintf(namefile, "./OUTPUT/BOX%d/STRESSES/angmom.dat",k);
        fid = fopen(namefile,"a");      //Append to end of file

        //Calculate the total angular momentum

        double Lx = 0.;
        double Ly = 0.;
        double Lz = 0.;
	double A[3], B[3], C[3];
	int ns = box[k].boxns;

	//Get the net linear momentum to subtract out when computing angular momentum
  double pvx = 0.0;
    double pvy = 0.0;
    double pvz = 0.0;
          double xc = 0.0;
          double yc = 0.0;
          double zc = 0.0;
          double mc = 0.0;

  /* ----------------------------------------------- */
  /* Accumulate the momentum and COM in each         */
  /* direaction.                                     */
  /* ----------------------------------------------- */
    for(int i=0; i<box[k].boxns; i++) {
      double mass = pott[k][i].mas;
      pvx += vv[k][i].x * mass;
      pvy += vv[k][i].y * mass;
      pvz += vv[k][i].z * mass;
                  xc += atnopbc[k][i].x * mass;
                  yc += atnopbc[k][i].y * mass;
                  zc += atnopbc[k][i].z * mass;
                  mc += pott[k][i].mas;
    }

  /* ----------------------------------------------- */
  /* Determine the average net momentum in each      */
  /* direction.                                      */
  /* ----------------------------------------------- */
    pvx /= (1.0*box[k].boxns);
    pvy /= (1.0*box[k].boxns);
    pvz /= (1.0*box[k].boxns);

	double sumx=0;

        for(int i=0; i<ns; i++)
        {
                double mass = pott[k][i].mas;
                double rx   = atom[k][i].x;
                double ry   = atom[k][i].y;
                double rz   = atom[k][i].z;
                //Calculate the total angular momentum
                A[0] = rx;
                A[1] = ry;
                A[2] = rz;
                B[0] = mass * vv[k][i].x - pvx;	//Subtract out the net linear mom.
                B[1] = mass * vv[k][i].y - pvy;
                B[2] = mass * vv[k][i].z - pvz;
		sumx+=B[0];
                CROSS(A,B,C);
                Lx += C[0];
                Ly += C[1];
                Lz += C[2];
        }
	
        //Calculate the simulation time.
        double tt;
        if (flag==0)
                tt = sim.dt * 1.0*icycle;
        else if (flag==1 && (sim.ID2==3 || sim.ID2==5))
                tt= sim.dt*sim.cyc_eq+ sim.dtlong*1.0*icycle; //multiple time step
        else
                tt = sim.dt * 1.0*(sim.cyc_eq+icycle);


        double angmom = sqrt(Lx*Lx+Ly*Ly+Lz*Lz);

        fprintf(fid, "%12.7f\t%12.7f\t\n", tt, angmom);
        fclose(fid);


	//Now compute and output the magnitude of the net linear momentum.
        sprintf(namefile2, "./OUTPUT/BOX%d/STRESSES/linmom.dat",k);
        fid2 = fopen(namefile2,"a");      //Append to end of file

        double linmom = sqrt(pvx*pvx+pvy*pvy+pvz*pvz);

	fprintf(fid2, "%12.7f\t%12.7f\t\n", tt, linmom);
	fclose(fid2);
}

#endif //ZEROAM

#endif //ZHOU
