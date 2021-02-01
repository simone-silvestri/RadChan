
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

using namespace std;


/****** CALCULATING THE NARROW-BANDS CORRELATED-K ******/

void correlatedK(deque<PartSpecBand> &bands, int totbands)
{
	double pwr = 0.1;
	int iaddo;
	double kmax = -1000,kmin = 1e10;
	double ktmp, pwrk[n_pwrk+1],ki,kpwri;
	double kBarLbl, kBarCk;
	double maxError = 0;
	double meanError = 0;
	int bandError[2];

	// find minumum and maximum kappa values in narrow-band

	for (int i=0; i<bands.size(); i++)
	{
		for (int j=0; j<bands[i].numbands; j++)
		{
			for (int h=0; h<bands[i].narrband[j].kvn.size(); h++)
			{
				ktmp = bands[i].narrband[j].kvn[h];
				if (ktmp > kmax)
					kmax = ktmp;
				if (ktmp < kmin)
					kmin = ktmp;
			}
			bands[i].narrband[j].kmin = kmin;
			bands[i].narrband[j].kmax = kmax;
			kmax= - 1000;
			kmin= 1e10;
		}
	}

	/***********  START WITH K-DISTRIBUTION **********/

	// loop on every narrow-band in every wide-band to redistribute k
	for (int i=0; i<bands.size(); i++)
	{
		for (int j=0; j<bands[i].numbands; j++)
		{
			double pwrk_min = pow( bands[i].narrband[j].kmin, pwr);
			double pwrk_max = pow( bands[i].narrband[j].kmax, pwr);

			double pwrk_step = (pwrk_max - pwrk_min) / (n_pwrk - 1);

			// create k-distribution with logarithmic spacing
			for ( int h=0; h<n_pwrk+1; h++)
			{
				pwrk[h] = (h)*pwrk_step + pwrk_min;
				bands[i].narrband[j].kk[h] = pow( pwrk[h], 1/pwr );
			}

			// initialize probability distribution to 0
			for (int k=0; k<n_pwrk+1; k++)
			{
				bands[i].narrband[j].ff[k] = 0.0;
				bands[i].narrband[j].gg[k] = 0.0;
			}

			// start redistributing now
			for ( int h=1; h<bands[i].narrband[j].kvn.size(); h++)
			{
				ki = (bands[i].narrband[j].kvn[h]+bands[i].narrband[j].kvn[h-1])/2;
				kpwri = pow( ki, pwr);
				iaddo = MAX(0 , (int) ((kpwri - pwrk_min)/pwrk_step+1));

				bands[i].narrband[j].ff[iaddo] = bands[i].narrband[j].ff[iaddo] + (bands[i].narrband[j].wvn[h]-bands[i].narrband[j].wvn[h-1])/
						(bands[i].narrband[j].wv_right-bands[i].narrband[j].wv_left);

			}// closed loop into the same narrow-band

			// last bin is empty
			int ng = n_pwrk - 1;

			//correct k-values by half a bin upwards
			for (int h=0; h<n_pwrk; h++)
				bands[i].narrband[j].kk[h] = 0.5*(bands[i].narrband[j].kk[h]+bands[i].narrband[j].kk[h+1]);

			//calculate raw g-function
			bands[i].narrband[j].gg[ng] = 1.0;
			for (int h=ng; h>0; h--)
				bands[i].narrband[j].gg[h-1] = bands[i].narrband[j].gg[h] - bands[i].narrband[j].ff[h];

			/****** Now calculate the exactness of the c-k distribution **************/

			kBarLbl = 0.0;
			kBarCk  = 0.0;
			for (int h=0; h<bands[i].narrband[j].wvn.size() ; h++)
				kBarLbl = kBarLbl + bands[i].narrband[j].kvn[h]/(double) bands[i].narrband[j].wvn.size();
			for (int h=0; h<n_pwrk ; h++)
				kBarCk = kBarCk + bands[i].narrband[j].kk[h]*(bands[i].narrband[j].ff[h]);;

			maxError = MAX(maxError, (kBarLbl-kBarCk)/kBarLbl*100);
			meanError = meanError + (kBarLbl-kBarCk)/kBarLbl*100/totbands;
			if (maxError == (kBarLbl-kBarCk)/kBarLbl*100)
			{
				bandError[1] = i;
				bandError[2] = j;
			}

			// calculate k on quadrature point gq

			int t = 0;
			for (int r=0; r<nQuad; r++)
			{
				RESTART:
				if(bands[i].narrband[j].gg[t+1] >= gq[r])
				{
					bands[i].narrband[j].kq[r] = bands[i].narrband[j].kk[t]+(bands[i].narrband[j].kk[t+1]-bands[i].narrband[j].kk[t])*
							(gq[r] - bands[i].narrband[j].gg[t])/(bands[i].narrband[j].gg[t+1] - bands[i].narrband[j].gg[t]);
				}
				else
				{
					t=t+1;
					goto RESTART;
				}
			}

		}//close loop on narrow-bands
	}//closed loop on widebands

	printf("Max error in -> Wide-band %d, Narrow-band %d, max error %lf %%, mean error %lf %%\n",bandError[1],bandError[2],maxError,meanError);

}

void correlatedKT(TempBand *tempBand)
{
	double pwr = 0.1;
	int iaddo;
	double kmax = -1000,kmin = 1e10;
	double ktmp, pwrk[n_pwrk+1],ki,kpwri;
	double kBarLbl, kBarCk, kBarCkq;
	double tBarLbl, tBarCk, tBarCkq;
	double aBarLbl, aBarCk, aBarCkq;
	double taBarLbl, taBarCk, taBarCkq;

	// find minumum and maximum kappa values in narrow-band

	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			for (int h=0; h<tempBand[i].narrband[j].size; h++)
			{
				ktmp = tempBand[i].narrband[j].kvn[h];
				if (ktmp > kmax)
					kmax = ktmp;
				if (ktmp < kmin)
					kmin = ktmp;
			}
			tempBand[i].narrband[j].kmin = kmin;
			tempBand[i].narrband[j].kmax = kmax;
			kmax= - 1000;
			kmin= 1e10;
		}
	}

	/***********  START WITH K-DISTRIBUTION **********/

	// loop on every narrow-band in every temperature to redistribute k
	for (int i=0; i<nT; i++)
	{
		sprintf(tempBand[i].file_err,"errors/error%d.txt",(int) tempBand[i].T);
		sprintf(tempBand[i].file_par,"errors/param%d.txt",(int) tempBand[i].T);

		FILE *fp = fopen(tempBand[i].file_err,"wt");
		FILE *fp1= fopen(tempBand[i].file_par,"wt");
		printf("Re-distributing temperature %lf\n",tempBand[i].T);
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			double pwrk_min = pow( tempBand[i].narrband[j].kmin, pwr);
			double pwrk_max = pow( tempBand[i].narrband[j].kmax, pwr);

			double pwrk_step = (pwrk_max - pwrk_min) / (n_pwrk - 1);

			// create k-distribution with logarithmic spacing
			for ( int h=0; h<n_pwrk+1; h++)
			{
				pwrk[h] = (h)*pwrk_step + pwrk_min;
				tempBand[i].narrband[j].kk[h] = pow( pwrk[h], 1/pwr );
			}

			// initialize probability distribution to 0
			for (int k=0; k<n_pwrk+1; k++)
			{
				tempBand[i].narrband[j].ff[k] = 0.0;
				tempBand[i].narrband[j].gg[k] = 0.0;
			}

			// start redistributing now
			for ( int h=1; h<tempBand[i].narrband[j].size; h++)
			{
				ki = (tempBand[i].narrband[j].kvn[h]+tempBand[i].narrband[j].kvn[h-1])/2;
				kpwri = pow( ki, pwr);
				iaddo = MAX(0 , (int) ((kpwri - pwrk_min)/pwrk_step+1));

				tempBand[i].narrband[j].ff[iaddo] = tempBand[i].narrband[j].ff[iaddo] + (tempBand[i].narrband[j].wvn[h]-tempBand[i].narrband[j].wvn[h-1])/
						(tempBand[i].narrband[j].wv_right-tempBand[i].narrband[j].wv_left);

			}// closed loop into the same narrow-band

			// last bin is empty
			int ng = n_pwrk - 1;

			//correct k-values by half a bin upwards
			for (int h=0; h<n_pwrk; h++)
				tempBand[i].narrband[j].kk[h] = 0.5*(tempBand[i].narrband[j].kk[h]+tempBand[i].narrband[j].kk[h+1]);

			//calculate raw g-function
			tempBand[i].narrband[j].gg[ng] = 1.0;
			for (int h=ng; h>0; h--)
				tempBand[i].narrband[j].gg[h-1] = tempBand[i].narrband[j].gg[h] - tempBand[i].narrband[j].ff[h];


			// calculate k on quadrature point gq

			int t = 0;
			for (int r=0; r<nQuad; r++)
			{
				RESTART:
				if(tempBand[i].narrband[j].gg[t+1] >= gq[r])
				{
					tempBand[i].narrband[j].kq[r] = tempBand[i].narrband[j].kk[t]+(tempBand[i].narrband[j].kk[t+1]-tempBand[i].narrband[j].kk[t])*
							(gq[r] - tempBand[i].narrband[j].gg[t])/(tempBand[i].narrband[j].gg[t+1] - tempBand[i].narrband[j].gg[t]);
					if(tempBand[i].narrband[j].gg[t+1]==tempBand[i].narrband[j].gg[t])
						tempBand[i].narrband[j].kq[r] = tempBand[i].narrband[j].kk[t+1];
				}
				else
				{
					t=t+1;
					goto RESTART;
				}
			}

			/****** Now calculate the exactness of the c-k distribution **************/

			kBarLbl = 0.0;
			kBarCk  = 0.0;
			kBarCkq = 0.0;
			tBarLbl = 0.0;
			tBarCk = 0.0;
			tBarCkq = 0.0;
			aBarLbl = 0.0;
			aBarCk = 0.0;
			aBarCkq = 0.0;
			taBarLbl = 0.0;
			taBarCk = 0.0;
			taBarCkq = 0.0;

			for (int h=1; h<tempBand[i].narrband[j].size ; h++)
			{
				kBarLbl = kBarLbl + tempBand[i].narrband[j].kvn[h] * dwv / ( tempBand[i].narrband[j].wv_right - tempBand[i].narrband[j].wv_left );//(double) tempBand[i].narrband[j].size;
				tBarLbl = tBarLbl + exp(-tempBand[i].narrband[j].kvn[h]*0.5) * dwv / ( tempBand[i].narrband[j].wv_right - tempBand[i].narrband[j].wv_left ); //(double) tempBand[i].narrband[j].size;
				aBarLbl = aBarLbl + (1 - exp(-(tempBand[i].narrband[j].kvn[h]+tempBand[i].narrband[j].kvn[h-1])/2*0.5) )* dwv / ( tempBand[i].narrband[j].wv_right - tempBand[i].narrband[j].wv_left ); //(double) tempBand[i].narrband[j].size;
				taBarLbl = taBarLbl + exp(-tempBand[i].narrband[j].kvn[h]) * (1 - exp(-tempBand[i].narrband[j].kvn[h]) )* dwv / ( tempBand[i].narrband[j].wv_right - tempBand[i].narrband[j].wv_left ); //(double) tempBand[i].narrband[j].size;
			}
			for (int h=1; h<n_pwrk ; h++)
			{
				kBarCk = kBarCk + (tempBand[i].narrband[j].kk[h])*(tempBand[i].narrband[j].ff[h]);
				tBarCk = tBarCk + exp(-tempBand[i].narrband[j].kk[h]*0.5)*(tempBand[i].narrband[j].ff[h]);
				aBarCk = aBarCk + (1 -exp(-tempBand[i].narrband[j].kk[h]))*(tempBand[i].narrband[j].ff[h]);
				taBarCk = taBarCk +  exp(-tempBand[i].narrband[j].kk[h])*(1 -exp(-tempBand[i].narrband[j].kk[h]))*(tempBand[i].narrband[j].ff[h]);
			}
			for (int h=0; h<nQuad ; h++)
			{
				kBarCkq = kBarCkq + tempBand[i].narrband[j].kq[h]*(wq[h]);
				tBarCkq = tBarCkq + exp(-tempBand[i].narrband[j].kq[h]*0.5)*(wq[h]);
				aBarCkq = aBarCkq + (1 -exp(-tempBand[i].narrband[j].kq[h]))*(wq[h]);
				taBarCkq = taBarCkq + exp(-tempBand[i].narrband[j].kq[h])*(1 -exp(-tempBand[i].narrband[j].kq[h]))*(wq[h]);
			}

			fprintf(fp,"%d %le %le %le %le %le %le %le %le\n",j,(kBarLbl-kBarCk)/kBarLbl*100,(kBarLbl-kBarCkq)/kBarLbl*100,
					(tBarLbl-tBarCk)/tBarLbl*100,(tBarLbl-tBarCkq)/tBarLbl*100,(aBarLbl-aBarCk)/aBarLbl*100,(aBarLbl-aBarCkq)/aBarLbl*100,
					(taBarLbl-taBarCk)/taBarLbl*100,(taBarLbl-taBarCkq)/taBarLbl*100);
			fprintf(fp1,"%d %le %le %le %le %le %le %le %le %le %le %le %le\n",j,kBarLbl,kBarCk,kBarCkq,tBarLbl,tBarCk,tBarCkq,aBarLbl,aBarCk,aBarCkq
					,taBarLbl,taBarCk,taBarCkq);
		}//close loop on narrow-bands
		fclose(fp);
		fclose(fp1);
	}//closed loop on Temperature

}

void correlatedKP(PartSpec2 *particles)
{
	double pwr = 0.1;
	int iaddok, iaddos;
	double kmax = -1000,kmin = 1e10;
	double smax = -1000,smin = 1e10;
	double ktmp, kpwrk[n_pwrk+1],ki,kpwri;
	double stmp, spwrk[n_pwrk+1],si,spwri;
	double kBarLbl, kBarCk, kBarCkq;
	double sBarLbl, sBarCk, sBarCkq;

	// find minumum and maximum kappa values in narrow-band

	for (int j = 0; j < particles->totbands; j++)
	{
		for (int h = 0; h < particles->narr_abs[j].size; h++)
		{
			ktmp = particles->narr_abs[j].kvn[h];
			stmp = particles->narr_sca[j].kvn[h];
			if (ktmp > kmax)
				kmax = ktmp;
			if (ktmp < kmin)
				kmin = ktmp;
			if (stmp > smax)
				smax = stmp;
			if (stmp < smin)
				smin = stmp;
		}
		particles->narr_abs[j].kmin = kmin;
		particles->narr_abs[j].kmax = kmax;
		particles->narr_sca[j].kmin = smin;
		particles->narr_sca[j].kmin = smax;
		kmax= - 1000;
		kmin= 1e10;
	}

	/***********  START WITH K-DISTRIBUTION **********/


	sprintf(particles->file_err,"errors/error-part.txt");
	sprintf(particles->file_par,"errors/param-part.txt");

	FILE *fp = fopen(particles->file_err,"wt");
	FILE *fp1= fopen(particles->file_par,"wt");
	printf("Re-distributing particles from %s\n",particles->file);
	for (int j=0; j < particles->totbands; j++)
	{
		double kpwrk_min = pow(particles->narr_abs[j].kmin, pwr);
		double kpwrk_max = pow(particles->narr_abs[j].kmax, pwr);
		double spwrk_min = pow(particles->narr_sca[j].kmin, pwr);
		double spwrk_max = pow(particles->narr_sca[j].kmax, pwr);

		double kpwrk_step = (kpwrk_max - kpwrk_min) / (n_pwrk - 1);
		double spwrk_step = (spwrk_max - spwrk_min) / (n_pwrk - 1);

		// create k-distribution with logarithmic spacing
		for ( int h=0; h<n_pwrk+1; h++)
		{
			kpwrk[h] = (h)*kpwrk_step + kpwrk_min;
			spwrk[h] = (h)*spwrk_step + spwrk_min;

			particles->narr_abs[j].kk[h] = pow( kpwrk[h], 1/pwr );
			particles->narr_sca[j].kk[h] = pow( spwrk[h], 1/pwr );

		}

		// initialize probability distribution to 0
		for (int k=0; k<n_pwrk+1; k++)
		{
			particles->narr_abs[j].ff[k] = 0.0;
			particles->narr_abs[j].gg[k] = 0.0;
			particles->narr_sca[j].ff[k] = 0.0;
			particles->narr_sca[j].gg[k] = 0.0;
		}

		// start redistributing now
		for ( int h=1; h<particles->narr_abs[j].size; h++)
		{
			ki = (particles->narr_abs[j].kvn[h]+particles->narr_abs[j].kvn[h-1])/2;
			kpwri = pow( ki, pwr);
			iaddok = MAX(0 , (int) ((kpwri - kpwrk_min)/kpwrk_step+1));

			particles->narr_abs[j].ff[iaddok] = particles->narr_abs[j].ff[iaddok] + (particles->narr_abs[j].wvn[h]-particles->narr_abs[j].wvn[h-1])/
					(particles->narr_abs[j].wv_right-particles->narr_abs[j].wv_left);

			si = (particles->narr_sca[j].kvn[h]+particles->narr_sca[j].kvn[h-1])/2;
			spwri = pow( si, pwr);
			iaddos = MAX(0 , (int) ((spwri - spwrk_min)/spwrk_step+1));

			particles->narr_sca[j].ff[iaddos] = particles->narr_sca[j].ff[iaddos] + (particles->narr_sca[j].wvn[h]-particles->narr_sca[j].wvn[h-1])/
					(particles->narr_sca[j].wv_right-particles->narr_sca[j].wv_left);

		}// closed loop into the same narrow-band

		// last bin is empty
		int ng = n_pwrk - 1;

		//correct k-values by half a bin upwards
		for (int h=0; h<n_pwrk; h++) {
			particles->narr_abs[j].kk[h] = 0.5*(particles->narr_abs[j].kk[h]+particles->narr_abs[j].kk[h+1]);
			particles->narr_sca[j].kk[h] = 0.5*(particles->narr_sca[j].kk[h]+particles->narr_sca[j].kk[h+1]);
		}
		//calculate raw g-function
		particles->narr_abs[j].gg[ng] = 1.0;
		particles->narr_sca[j].gg[ng] = 1.0;
		for (int h=ng; h>0; h--) {
			particles->narr_abs[j].gg[h-1] = particles->narr_abs[j].gg[h] - particles->narr_abs[j].ff[h];
			particles->narr_sca[j].gg[h-1] = particles->narr_sca[j].gg[h] - particles->narr_sca[j].ff[h];
		}

		// calculate k on quadrature point gq

		int t = 0;
		for (int r=0; r<nQuad; r++)
		{
			RESTARTK:
			if(particles->narr_abs[j].gg[t+1] >= gq[r])
			{
				particles->narr_abs[j].kq[r] = particles->narr_abs[j].kk[t]+(particles->narr_abs[j].kk[t+1]-particles->narr_abs[j].kk[t])*
						(gq[r] - particles->narr_abs[j].gg[t])/(particles->narr_abs[j].gg[t+1] - particles->narr_abs[j].gg[t]);
				if(particles->narr_abs[j].gg[t+1]==particles->narr_abs[j].gg[t])
					particles->narr_abs[j].kq[r] = particles->narr_abs[j].kk[t+1];
			}
			else
			{
				t=t+1;
				goto RESTARTK;
			}
		}
		// calculate s on quadrature point gq
		t = 0;
		for (int r=0; r<nQuad; r++)
		{
			RESTARTS:
			if(particles->narr_sca[j].gg[t+1] >= gq[r])
			{
				particles->narr_sca[j].kq[r] = particles->narr_sca[j].kk[t]+(particles->narr_sca[j].kk[t+1]-particles->narr_sca[j].kk[t])*
						(gq[r] - particles->narr_sca[j].gg[t])/(particles->narr_sca[j].gg[t+1] - particles->narr_sca[j].gg[t]);
				if(particles->narr_sca[j].gg[t+1]==particles->narr_sca[j].gg[t])
					particles->narr_sca[j].kq[r] = particles->narr_sca[j].kk[t+1];
			}
			else
			{
				t=t+1;
				goto RESTARTS;
			}
		}
		/****** Now calculate the exactness of the c-k distribution **************/

		kBarLbl = 0.0;
		kBarCk  = 0.0;
		kBarCkq = 0.0;

		sBarLbl = 0.0;
		sBarCk  = 0.0;
		sBarCkq = 0.0;

		for (int h=1; h<particles->narr_abs[j].size ; h++)
		{
			kBarLbl += particles->narr_abs[j].kvn[h] * dwv / ( particles->narr_abs[j].wv_right - particles->narr_abs[j].wv_left );//(double) tempBand[i].narrband[j].size;
			sBarLbl += particles->narr_sca[j].kvn[h] * dwv / ( particles->narr_sca[j].wv_right - particles->narr_sca[j].wv_left );//(double) tempBand[i].narrband[j].size;
		}
		for (int h=1; h<n_pwrk ; h++)
		{
			kBarCk += (particles->narr_abs[j].kk[h])*(particles->narr_abs[j].ff[h]);
			sBarCk += (particles->narr_sca[j].kk[h])*(particles->narr_sca[j].ff[h]);
		}
		for (int h=0; h<nQuad ; h++)
		{
			kBarCkq += particles->narr_abs[j].kq[h]*(wq[h]);
			sBarCkq += particles->narr_sca[j].kq[h]*(wq[h]);
		}

		fprintf(fp,"%le %le %le %le %le\n",particles->narr_abs[j].wvc,(kBarLbl-kBarCk)/kBarLbl*100,(kBarLbl-kBarCkq)/kBarLbl*100,(sBarLbl-sBarCk)/sBarLbl*100,(sBarLbl-sBarCkq)/sBarLbl*100);
		fprintf(fp1,"%le %le %le %le %le %le %le\n",particles->narr_abs[j].wvc,kBarLbl,kBarCk,kBarCkq,sBarLbl,sBarCk,sBarCkq);
	}//close loop on narrow-bands
	fclose(fp);
	fclose(fp1);

}
