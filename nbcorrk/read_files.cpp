
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
#include <deque>
#include <iostream>

using namespace std;
double lininterp(double y1, double y2, double x1, double x2, double x);
/****** READING AND CALCULATING THE WIDE-BANDS ******/

double absWideBand(deque<PartSpecBand> &bands, char *namefile, double T)
{
	double wv,kv,Ibv;

	FILE *fp = fopen(namefile,"r");

	int i = 0, n = 0;
	double kvprev = 0.0,wvprev = 0.0,total = 0.0, wvcu=wvtoll+1.;
	while (fscanf(fp,"%lf %lf \n",&wv,&kv) != EOF)
	{

//      kv = kv*H+0.02;
		if (kv > toll)
		{
			if (kvprev <= toll && wvcu >= wvtoll)
			{
				n++;
				PartSpecBand band;
				band.kv.push_back(kv*100);  //remember to scan it as m**(-1) and not cm**(-1)
				band.wv.push_back(wv);
				bands.push_back(band);
			}
			else
			{
				bands.back().kv.push_back(kv*100);  //remember to scan it as m**(-1) and not cm**(-1)
				bands.back().wv.push_back(wv);
			}
			wvcu = 0.0;
		}
		else
		{
			if (wvcu < wvtoll)
			{
				bands.back().kv.push_back(kv*100);  //remember to scan it as m**(-1) and not cm**(-1)
				bands.back().wv.push_back(wv);
				wvcu = wvcu + dwv;
			}
		}

		Ibv = I_black( T, (wv + wvprev)/2 );
		total = total + (kv + kvprev)/2 * Ibv * (wv - wvprev);

		kvprev = kv;
		wvprev = wv;
		i++;

	};

	fclose(fp);
	printf("Total intensity : %lf, Total Wide bands : %d \n",total,n);
	return total;

}

void readFillT(TempBand *tempBand, deque<PartSpecBand> &bands)
{
	double wv,kv,dummy;

	for (int i=0; i<nT; i++)
	{
		FILE *fp = fopen(tempBand[i].file,"r");

		printf("Reading data from temperature: %lf\n",tempBand[i].T);
		int k = 0;
		int t = 0;
		int flag = 0;
		tempBand[i].kPr = 0;
		while (fscanf(fp,"%lf %lf \n",&wv,&kv) != EOF)
		{

//         kv = kv*H+0.02;
			//remember to calculate kP in m**(-1), therefore also wvn must be changed in m**(-1)
			tempBand[i].kPr = tempBand[i].kPr + pi * kv *100 * I_black(tempBand[i].T,wv) * dwv * 100 / ( stefan * pow(tempBand[i].T, 4));

			if (wv == tempBand[i].narrband[k].wvn[t])
			{
				tempBand[i].narrband[k].kvn[t] = kv*100; //remember to scan it as m**(-1) and not cm**(-1)
				t = t+1;
				if ( t==tempBand[i].narrband[k].size )
				{
					t = 0;
					k = k+1;
					if ( k == tempBand[i].totbands )
						break;
				}
			}

		};
		fclose(fp);

	}

}

void readFillN(TempBand *tempBand, deque<PartSpecBand> &bands)
{
	double wv,kv,dummy;

	for (int i=0; i<nT; i++)
	{
		FILE *fp = fopen(tempBand[i].file,"r");

		printf("Reading data from temperature: %lf\n",tempBand[i].T);
		int k = 0;
		int t = 0;
		int flag = 0;
		tempBand[i].kPr = 0;
		while (fscanf(fp,"%lf %lf \n",&wv,&kv) != EOF)
		{

         kv *= H;
			//remember to calculate kP in m**(-1), therefore also wvn must be changed in m**(-1)
			tempBand[i].kPr = tempBand[i].kPr + pi * kv *100 * I_black(tempBand[i].T,wv) * dwv * 100 / ( stefan * pow(tempBand[i].T, 4));

			if (wv == tempBand[i].narrband[k].wvn[t])
			{
				tempBand[i].narrband[k].kvn[t] = kv*100; //remember to scan it as m**(-1) and not cm**(-1)
				t = t+1;
				if ( t==tempBand[i].narrband[k].size )
				{
					t = 0;
					k = k+1;
					if ( k == tempBand[i].totbands )
						break;
				}
			}

		};
		fclose(fp);

	}

}
void readFillP(PartSpec2 *particles, deque<PartSpecBand> &bands)
{
	double wv,ev,sv,dummy1,dummy2;
	int dummy;
	int idx;
	PartSpec partspec;
	FILE *fp = fopen(particles->file,"r");
	printf("Reading particle data from file: %s\n",particles->file);

	while (fscanf(fp,"%d %le %le %le %le %le\n",&dummy,&wv,&sv,&ev,&dummy1,&dummy2) != EOF)
	{
		partspec.wv.push_back(wv);
		partspec.kv.push_back(ev-sv);
		partspec.sv.push_back(sv);
	}
	fclose(fp);
	for(int i = 0; i < particles->totbands; i++)
	{
		for(int j = 0; j < particles->narr_abs[i].size; j++)
		{
			// find index of partspec where wv lies and linear interpolation
			if(particles->narr_abs[i].wvn[j] < partspec.wv[0])
			{
				particles->narr_abs[i].kvn[j] = 0;
				particles->narr_sca[i].kvn[j] = 0;
			}
			else if(particles->narr_abs[i].wvn[j] > partspec.wv[partspec.wv.size()-1])
			{
				particles->narr_abs[i].kvn[j] = partspec.kv[partspec.wv.size()-1];
				particles->narr_sca[i].kvn[j] = partspec.sv[partspec.wv.size()-1];
			}
			else
			{
				for(int k = 1; k < partspec.wv.size(); k++)
				{
					if((particles->narr_abs[i].wvn[j] < partspec.wv[k]) && (particles->narr_abs[i].wvn[j] > partspec.wv[k-1])) {
						idx = k;
						break;
					}
				}
				particles->narr_abs[i].kvn[j] = lininterp(partspec.kv[idx-1],partspec.kv[idx],partspec.wv[idx-1],partspec.wv[idx],particles->narr_sca[i].wvn[j]);
				particles->narr_sca[i].kvn[j] = lininterp(partspec.sv[idx-1],partspec.sv[idx],partspec.wv[idx-1],partspec.wv[idx],particles->narr_sca[i].wvn[j]);
			}
		}
	}

}

void readFillH(Phase *phi, deque<PartSpecBand> &bands)
{
	double wv,tv,phiv;
	double tv1 = 0;
	int idx;
	int t = 0;
	int m = 0;
	Phaset phit;
	FILE *fp = fopen(phi->file,"r");
	printf("Reading phase data from file: %s\n",phi->file);
	while ((fscanf(fp,"%le %le %le\n",&wv,&tv,&phiv) != EOF))
	{
		if(t==180) {
			phit.wv.push_back(wv);
			t = 0;
		}
		phit.phiv[t].push_back(phiv);
		t += 1;
	}
	fclose(fp);
	for(int t = 0; t < nA; t++)
	{
		phi->t[t] = t*pi/(nA*1.0);
		for(int i = 0; i < phi->totbands; i++)
		{
			// find index of partspec where wv lies and linear interpolation
			if(phi->wvc[i] < phit.wv[0])
			{
				phi->p[t][i] = phit.phiv[t][0];
			}
			else if(phi->wvc[i] > phit.wv[phit.wv.size()-1])
			{
				phi->p[t][i] = phit.phiv[t][phit.wv.size()-1];
			}
			else
			{
				for(int k = 1; k < phit.wv.size(); k++)
				{
					if((phi->wvc[i] < phit.wv[k]) && (phi->wvc[i] > phit.wv[k-1])) {
						idx = k;
						break;
					}
				}
				phi->p[t][i] = lininterp(phit.phiv[t][idx-1],phit.phiv[t][idx],phit.wv[idx-1],phit.wv[idx],phi->wvc[i]);
			}
		}
	}
}

double lininterp(double y1, double y2, double x1, double x2, double x)
{
	return (y2 - y1) / (x2 - x1) * (x - x1) + y1;
}
