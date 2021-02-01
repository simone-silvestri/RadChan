
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


void allocT(TempBand *tempBand, deque<PartSpecBand> &bands, int totbands)
{
	for (int i=0; i<nT; i++)
	{
		tempBand[i].T = Tinit + DT*i;
		tempBand[i].allocNarrBands(totbands);
		tempBand[i].totbands = totbands;
		int k=0;
		for (int j=0; j<bands.size(); j++)
		{
			for (int h=0; h<bands[j].numbands; h++)
			{
				tempBand[i].narrband[k].allocVar(bands[j].narrband[h].wvn.size());
				tempBand[i].narrband[k].size = bands[j].narrband[h].wvn.size();
				tempBand[i].narrband[k].wv_left  = bands[j].narrband[h].wv_left;
				tempBand[i].narrband[k].wv_right = bands[j].narrband[h].wv_right;
				tempBand[i].narrband[k].wvc = bands[j].narrband[h].wvc;

				for (int t=0; t<bands[j].narrband[h].wvn.size(); t++)
				{
					tempBand[i].narrband[k].wvn[t] = bands[j].narrband[h].wvn[t];
				}

				k=k+1;
			}
		}
		if(tempBand[i].T >= 1000)
		{
			sprintf(tempBand[i].file,"CO2-data/CO2_%dK_HITRAN.txt",(int) tempBand[i].T);
		}
		else
		{
			sprintf(tempBand[i].file,"CO2-data/CO2_%dK_HITRAN.txt",(int) tempBand[i].T);
		}
	}
}

void allocP(PartSpec2 *particles, deque<PartSpecBand> &bands, int totbands)
{
	particles->allocNarrBands(totbands);
	particles->totbands = totbands;
	int k = 0;
	for (int j=0; j<bands.size(); j++)
	{
		for (int h=0; h<bands[j].numbands; h++)
		{
			particles->narr_abs[k].allocVar(bands[j].narrband[h].wvn.size());
			particles->narr_abs[k].size = bands[j].narrband[h].wvn.size();
			particles->narr_abs[k].wv_left  = bands[j].narrband[h].wv_left;
			particles->narr_abs[k].wv_right = bands[j].narrband[h].wv_right;
			particles->narr_abs[k].wvc = bands[j].narrband[h].wvc;
			particles->narr_sca[k].allocVar(bands[j].narrband[h].wvn.size());
			particles->narr_sca[k].size = bands[j].narrband[h].wvn.size();
			particles->narr_sca[k].wv_left  = bands[j].narrband[h].wv_left;
			particles->narr_sca[k].wv_right = bands[j].narrband[h].wv_right;
			particles->narr_sca[k].wvc = bands[j].narrband[h].wvc;

			for (int t=0; t<bands[j].narrband[h].wvn.size(); t++)
			{
				particles->narr_abs[k].wvn[t] = bands[j].narrband[h].wvn[t];
				particles->narr_sca[k].wvn[t] = bands[j].narrband[h].wvn[t];
			}
			k=k+1;
		}
	}
	sprintf(particles->file,"particle-data/carbon.txt");
}

void allocN(PartSpec2 *particles, deque<PartSpecBand> &bands, int totbands)
{
	particles->allocNarrBands(totbands);
	particles->totbands = totbands;
	int k = 0;
	for (int j=0; j<bands.size(); j++)
	{
		for (int h=0; h<bands[j].numbands; h++)
		{
			particles->narr_abs[k].allocVar(bands[j].narrband[h].wvn.size());
			particles->narr_abs[k].size = bands[j].narrband[h].wvn.size();
			particles->narr_abs[k].wv_left  = bands[j].narrband[h].wv_left;
			particles->narr_abs[k].wv_right = bands[j].narrband[h].wv_right;
			particles->narr_abs[k].wvc = bands[j].narrband[h].wvc;
			particles->narr_sca[k].allocVar(bands[j].narrband[h].wvn.size());
			particles->narr_sca[k].size = bands[j].narrband[h].wvn.size();
			particles->narr_sca[k].wv_left  = bands[j].narrband[h].wv_left;
			particles->narr_sca[k].wv_right = bands[j].narrband[h].wv_right;
			particles->narr_sca[k].wvc = bands[j].narrband[h].wvc;

			for (int t=0; t<bands[j].narrband[h].wvn.size(); t++)
			{
				particles->narr_abs[k].wvn[t] = bands[j].narrband[h].wvn[t];
				particles->narr_sca[k].wvn[t] = bands[j].narrband[h].wvn[t];
			}
			k=k+1;
		}
	}
	sprintf(particles->file,"particle-data/carbon-nano.txt");
}

void allocH(Phase *phi, deque<PartSpecBand> &bands, int totbands)
{
	phi->totbands = totbands;
	phi->allocVar(totbands);
	int k = 0;
	for (int j=0; j<bands.size(); j++)
	{
		for (int h=0; h<bands[j].numbands; h++)
		{
			phi->wvc[k] = bands[j].narrband[h].wvc;
			k=k+1;
		}
	}
	sprintf(phi->file,"particle-data/carbon-scatt.txt");
}
