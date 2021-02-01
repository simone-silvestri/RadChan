
/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <stdio.h>
#include <stdlib.h>

#include "param.h"
#include "functions.h"

void calcProbT(TempBand *tempBand)
{
	for (int i=0; i<nT; i++)
	{
		printf("Calculating probability for temperature %lf \n",tempBand[i].T);
		tempBand[i].kP = 0;
		tempBand[i].kR = 0;
		tempBand[i].kM = 0;
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			// H2O wavenumber high Re

			// H2O wavenumber
			double T  = (tempBand[i].T - 802.66)/87.514; //For H2O
                        double coeff[11] = {0.013638, 0.10321, 0.11378, -0.67375, -1.2894, 1.1924, 2.6412, -0.22784, 1.4089, -1.18, 4.7927}; //For H2O
                        //double coeff[11] = {0.018195, 0.14863, 0.19381, -0.93129, -1.9522, 1.5128, 3.65762, -0.058086, 2.9417, -1.4145, 4.8413};
 
			// CO2 wavenumber
                 	//double T  = (tempBand[i].T - 788.37)/81.641; //For CO2 
                        //double coeff[11] = {-0.005967, -0.0027836, 0.078597, 0.045018, 0.20969, -0.19999, -0.85065, 0.33099, 4.6508, -0.89307, 4.4758};  //For CO2

         		// Particle wavenumber
         		//double T  = (tempBand[i].T - 810.21)/92.602; //For PART   
		        //double coeff[11] = {-0.017996,-0.073358,-0.0099447,-0.029945,-0.84971,-0.10266,2.801,1.7093,1.3732,-1.7954,3.7819}; //For PART


			double wh = coeff[0]*pow(T,10) + coeff[1]*pow(T,9) + coeff[2]*pow(T,8) + coeff[3]*pow(T,7) 
			          + coeff[4]*pow(T,6 ) + coeff[5]*pow(T,5) + coeff[6]*pow(T,4) + coeff[7]*pow(T,3) 
			          + coeff[8]*pow(T,2 ) + coeff[9]*pow(T,1) + coeff[10];


			// wavenumber for modeling	
			//double wh = 0.000001; //7.0*(3750.0/2900.0);

			tempBand[i].narrband[j].k_avg = 0;
			for (int h=1; h<n_pwrk; h++)
			{
				tempBand[i].kP += tempBand[i].narrband[j].kk[h]*I_black(tempBand[i].T,tempBand[i].narrband[j].wvc)*(tempBand[i].narrband[j].ff[h])*
						(tempBand[i].narrband[j].wv_right - tempBand[i].narrband[j].wv_left)*100*pi/(stefan * pow(tempBand[i].T,4));
				tempBand[i].narrband[j].k_avg += tempBand[i].narrband[j].kk[h]*(tempBand[i].narrband[j].ff[h]);
			
				if(tempBand[i].narrband[j].kk[h]>toll)	{
					tempBand[i].kR += 1.0/tempBand[i].narrband[j].kk[h]*I_deriv(tempBand[i].T,tempBand[i].narrband[j].wvc)*(tempBand[i].narrband[j].ff[h])*
								(tempBand[i].narrband[j].wv_right - tempBand[i].narrband[j].wv_left)*100*pi/(4*stefan*pow(tempBand[i].T,3));  }

				tempBand[i].kM += tempBand[i].narrband[j].kk[h]*I_black(tempBand[i].T,tempBand[i].narrband[j].wvc)*(tempBand[i].narrband[j].ff[h])*
						(tempBand[i].narrband[j].wv_right - tempBand[i].narrband[j].wv_left)*100*pi/(stefan * pow(tempBand[i].T,4)) * 
                                                  tempBand[i].narrband[j].kk[h]/wh * atan(wh/tempBand[i].narrband[j].kk[h]);

			}
		}
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			tempBand[i].narrband[j].cuprob = 0;
			tempBand[i].narrband[j].prob = 0;

			for (int h=1; h<n_pwrk; h++)
			{
				tempBand[i].narrband[j].prob += tempBand[i].narrband[j].kk[h]*I_black(tempBand[i].T,tempBand[i].narrband[j].wvc)*(tempBand[i].narrband[j].ff[h])*
						(tempBand[i].narrband[j].wv_right - tempBand[i].narrband[j].wv_left)*100*pi/(tempBand[i].kP * stefan * pow(tempBand[i].T,4));
			}
			if(j>0)
			{
				tempBand[i].narrband[j].cuprob = tempBand[i].narrband[j-1].cuprob+tempBand[i].narrband[j].prob;
			}
			else
			{
				tempBand[i].narrband[j].cuprob = tempBand[i].narrband[j].prob;
			}

		}

		// probability "new - 2"
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			for (int g=0; g<nQuad; g++)
			{
				tempBand[i].narrband[j].prob_new2[g] = 0;
				tempBand[i].narrband[j].cuprob_new2[g] = 0;
			}
			for (int g=0; g<nQuad; g++)
			{
				tempBand[i].narrband[j].prob_new2[g] = tempBand[i].narrband[j].kq[g] * I_black(tempBand[i].T,tempBand[i].narrband[j].wvc)*
						(tempBand[i].narrband[j].wv_right - tempBand[i].narrband[j].wv_left) * 100 * wq[g];
			}
			for (int g=0; g<nQuad; g++)
			{
				if( g ==0 )
				{
					tempBand[i].narrband[j].cuprob_new2[g] = tempBand[i].narrband[j].prob_new2[g];
				}
				else
				{
					tempBand[i].narrband[j].cuprob_new2[g] = tempBand[i].narrband[j].cuprob_new2[g-1] + tempBand[i].narrband[j].prob_new2[g];
				}
			}
			for (int g=0; g<nQuad; g++)
			{
				tempBand[i].narrband[j].cuprob_new2[g] /= tempBand[i].narrband[j].cuprob_new2[nQuad-1];
			}
		}

		//LBL probability
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			for (int h=0; h<tempBand[i].narrband[j].size; h++)
			{
				tempBand[i].narrband[j].prob_lbl[h] = 0;
				tempBand[i].narrband[j].cuprob_lbl[h] = 0;
			}
			for (int h=0; h<tempBand[i].narrband[j].size; h++)
			{
				tempBand[i].narrband[j].prob_lbl[h] = tempBand[i].narrband[j].kvn[h]*dwv*I_black(tempBand[i].T,tempBand[i].narrband[j].wvn[h])*100*pi/(tempBand[i].kPr * stefan * pow(tempBand[i].T,4));
			}
		}
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			for (int h=0; h<tempBand[i].narrband[j].size; h++)
			{
				if(j==0 && h==0)
				{
					tempBand[i].narrband[j].cuprob_lbl[h] = tempBand[i].narrband[j].prob_lbl[h];
				}
				else if(j!=0 && h==0)
				{
					tempBand[i].narrband[j].cuprob_lbl[h] = tempBand[i].narrband[j-1].cuprob_lbl[tempBand[i].narrband[j-1].size-1]+tempBand[i].narrband[j].prob_lbl[h];
				}
				else
				{
					tempBand[i].narrband[j].cuprob_lbl[h] = tempBand[i].narrband[j].cuprob_lbl[h-1]+tempBand[i].narrband[j].prob_lbl[h];
				}
			}
		}
	}
}

void calcProbP(PartSpec2 *particles)
{
	printf("Calculating probability for particles from %s \n",particles->file);

	particles->calcPlanckabs();
	for (int j = 0; j < particles->totbands; j++) {
		particles->narr_abs[j].k_avg = 0;
		particles->narr_sca[j].k_avg = 0;
		for (int h=1; h<n_pwrk; h++) {
			particles->narr_abs[j].k_avg += particles->narr_abs[j].kk[h]*(particles->narr_abs[j].ff[h]);
			particles->narr_sca[j].k_avg += particles->narr_sca[j].kk[h]*(particles->narr_sca[j].ff[h]);
		}
	}
	for(int i = 0; i < nT; i++)
	{
		for (int j = 0; j < particles->totbands; j++)
		{
			particles->narr_abs[j].cuprob[i] = 0;
			particles->narr_abs[j].prob[i] = 0;

			for (int h=1; h<n_pwrk; h++)
			{
				particles->narr_abs[j].prob[i] += particles->narr_abs[j].kk[h]*I_black(particles->T[i], particles->narr_abs[j].wvc)*(particles->narr_abs[j].ff[h])*
						(particles->narr_abs[j].wv_right - particles->narr_abs[j].wv_left)*100*pi/(particles->kP[i] * stefan * pow(particles->T[i],4));
			}
			if(j>0)
			{
				particles->narr_abs[j].cuprob[i] = particles->narr_abs[j-1].cuprob[i]+particles->narr_abs[j].prob[i];
			}
			else
			{
				particles->narr_abs[j].cuprob[i] = particles->narr_abs[j].prob[i];
			}

		}

		// probability "new - 2"
		for (int j = 0; j < particles->totbands; j++)
		{
			for (int g = 0; g < nQuad; g++)
			{
				particles->narr_abs[j].prob_new2[i][g] = 0;
				particles->narr_abs[j].cuprob_new2[i][g] = 0;
			}
			for (int g = 0; g < nQuad; g++)
			{
				particles->narr_abs[j].prob_new2[i][g] = particles->narr_abs[j].kq[g] * I_black(particles->T[i], particles->narr_abs[j].wvc)*
						(particles->narr_abs[j].wv_right - particles->narr_abs[j].wv_left) * 100 * wq[g];
			}
			for (int g = 0; g < nQuad; g++)
			{
				if( g == 0 )
				{
					particles->narr_abs[j].cuprob_new2[i][g] = particles->narr_abs[j].prob_new2[i][g];
				}
				else
				{
					particles->narr_abs[j].cuprob_new2[i][g] = particles->narr_abs[j].cuprob_new2[i][g-1] + particles->narr_abs[j].prob_new2[i][g];
				}
			}
			for (int g=0; g<nQuad; g++)
			{
				particles->narr_abs[j].cuprob_new2[i][g] /= particles->narr_abs[j].cuprob_new2[i][nQuad-1];
			}
		}
                //LBL probability
                for (int j = 0; j < particles->totbands; j++)
                {
                        for (int h=0; h < particles->narr_abs[j].size; h++)
                        {
                                particles->narr_abs[j].prob_lbl[h] = 0;
                                particles->narr_abs[j].cuprob_lbl[h] = 0;
                        }
                        for (int h=0; h < particles->narr_abs[j].size; h++)
                        {
                                particles->narr_abs[j].prob_lbl[h] = particles->narr_abs[j].kvn[h]*dwv*I_black(particles->T[i],particles->narr_abs[j].wvn[h])*100*pi/(particles->kP[i] * stefan * pow(particles->T[i],4));
                        }
                }
                for (int j=0; j<particles->totbands; j++)
                {
                        for (int h=0; h<particles->narr_abs[j].size; h++)
                        {
                                if(j==0 && h==0)
                                {
                                        particles->narr_abs[j].cuprob_lbl[h] = particles->narr_abs[j].prob_lbl[h];
                                }
                                else if(j!=0 && h==0)
                                {
                                        particles->narr_abs[j].cuprob_lbl[h] = particles->narr_abs[j-1].cuprob_lbl[particles->narr_abs[j-1].size-1]+particles->narr_abs[j].prob_lbl[h];
                                }
                                else
                                {
                                        particles->narr_abs[j].cuprob_lbl[h] = particles->narr_abs[j].cuprob_lbl[h-1]+particles->narr_abs[j].prob_lbl[h];
                                }
                        }
                }

	}
}

void calcProbH(Phase *phi)
{
	for(int i = 0; i < phi->totbands; i++)
	{
		phi->cuprob[0][i] = 0;
		for (int t = 1; t < nA; t++)
		{
			phi->prob[t][i] = 0.5 * (phi->t[t] - phi->t[t-1]) *
									(phi->p[t][i]+phi->p[t-1][i]) * 0.5 *
									(sin(phi->t[t]) + sin(phi->t[t-1])) * 0.5;
			phi->cuprob[t][i] = phi->cuprob[t-1][i] + phi->prob[t][i];
		}
		for (int t = 1; t < nA; t++)
		{
			phi->cuprob[t][i] /= phi->cuprob[nA-1][i];
		}
	}
}
