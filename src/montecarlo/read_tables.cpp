
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
#include "../param/param.h"
#include "definitions.h"
#include "functions.h"

myfloat I_blackC( myfloat T, myfloat nu);
extern myfloat prob[nT][nB], probg[nT][nB][nQ];

#if grey == 0
void readT(NarrowBand *narrBand, myfloat *Tnb, myfloat *kP, myfloat *prob_h, myfloat *probg_h, myfloat Tmax, myfloat *kappamax)
{
 
    myfloat dummy,dummy1;
    //read planck-absorption coefficient
    FILE *fp = fopen("properties/planck-mean.txt","r");
    for(int i=0; i<nT; i++)
    {
 	   fscanf(fp,readf" "reade" "reade" "reade"\n",&Tnb[i],&kP[i],&dummy,&dummy1);
           kP[i] *= chanHeight;
    }
    fclose(fp);


    //absorption coefficient at maximum temperature in the system
    int t = (int) ((Tmax - Tnb[0])/(Tnb[1]-Tnb[0]));
    *kappamax = (kP[t+1] - kP[t])/(Tnb[t+1]-Tnb[t])*
   	   	   	   (Tmax - Tnb[t])+kP[t];

    //read tables of narrow-bands and probability
    char *file[nB];
    for (int i=0; i<nB; i++)
    {
       file[i] = (char*)malloc(68*sizeof(char));
 	   sprintf(file[i],"properties/NarrBand%d.txt",i);
 	   fp = fopen(file[i],"r");
 	   fscanf(fp,reade"\t"reade"\t"reade"\n",&narrBand[i].wvl,&narrBand[i].wvr,&narrBand[i].wvc);
 	   for (int j=0; j<nT; j++)
 	   {
 		   fscanf(fp," "readf"\t",&dummy);
 		   fscanf(fp," "readf"\t",&dummy1);
		   fscanf(fp," "readf"\t",&narrBand[i].kavg);
 		   for(int g=0; g<nQ; g++)
 		   {
 			   fscanf(fp," "reade"\t",&narrBand[i].kq[j][g]);
                           narrBand[i].kq[j][g] *= chanHeight;
 		   }
 		   fscanf(fp,"\n");
 	   }
 	   fclose(fp);
 	   narrBand[i].idx = i;
    }

    // read probability and define the host dynamic


    int dummy3;
    fp = fopen("properties/prob_new2.txt","r");
    int n=0;
    int m=0;
    for (int i=0; i<nT; i++)
    {
 	   for (int j=0; j<nB; j++)
 	   {
 		   fscanf(fp,readf" \t %d \t "reade" \t",&dummy,&dummy3,&prob[i][j]);
 		   prob_h[n] = prob[i][j];
 		   n += 1;

 		   for (int g=0; g<nQ; g++)
 		   {
 			   fscanf(fp,reade"\t",&probg[i][j][g]);
 	 		   probg_h[m] = probg[i][j][g];
 	 		   m += 1;

 		   }
 		   fscanf(fp,"\n");
 	   }
    }
    fclose(fp);

    for (int i=0; i<nB; i++)
    {
	free(file[i]);
    }
}

#else
void readT(NarrowBand *narrBand, myfloat *Tnb, myfloat *kP, myfloat *prob_h, myfloat *probg_h, myfloat Tmax, myfloat *kappamax)
{

        //Planck absorption coefficient is the grey absorption coefficient
	for(int i=0; i<nT; i++)
	{
		kP[i] = abscoeff;
		Tnb[i] = Tmins*0.9 + (Tmax*1.1 - Tmins*0.9) * ((myfloat)i)/((myfloat)nT);
	}

	// absorption coefficient at maximum temperature in the system
	*kappamax = abscoeff;

	//read tables of narrow-bands and probability
	myfloat wvl[nB],wvr[nB];
	wvl[0] = 0;
	for (int i=1; i<nB; i++)
	{
		wvl[i] = wvl[i-1]+wvmax/nB;
	}
	wvr[0] = wvmax/nB;
	for (int i=1; i<nB; i++)
	{
		wvr[i] = wvr[i-1]+wvmax/nB;
	}
	for (int i=0; i<nB; i++)
	{
		narrBand[i].wvc = (wvl[i]+wvr[i])/2.;
		narrBand[i].wvl = wvl[i];
		narrBand[i].wvr = wvr[i];
		for (int j=0; j<nT; j++)
		{
			for(int g=0; g<nQ; g++)
			{
				narrBand[i].kq[j][g] = abscoeff;
			}
		}
		narrBand[i].idx = i;
	}

	// read probability and define the host dynamic

	myfloat Ib[nT][nB];

	for(int t = 0; t<nT; t++ )
	{
		for(int nb = 0; nb<nB; nb++ )
		{
			Ib[t][nb] = I_blackC(Tnb[t],narrBand[nb].wvc);
		}
	}
	int n=0;
	int m=0;
	for (int i=0; i<nT; i++)
	{
		prob[i][0] =  Ib[i][0]*(wvr[0]-wvl[0])*100;
		for (int j=1; j<nB; j++)
		{
			prob[i][j] = prob[i][j-1] + Ib[i][j]*(wvr[j]-wvl[j])*100;
		}
	}
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<nB; j++)
		{
			prob[i][j] = prob[i][j]/prob[i][nB-1];
			prob_h[n] = prob[i][j];
			n += 1;

			for (int g=0; g<nQ; g++)
			{
				probg[nT][nB][nQ] = g;
				probg_h[m] = probg[nT][nB][nQ];
				m += 1;
			}
		}
	}
}
#endif

/*  readS: function for calculating CDF of grey surfaces
 *  
 *  \param a
 *
 *
 *
 * 								            
 */ 								            


void readS(NarrowBand *narrBand, myfloat *Tnb, myfloat Tmax, myfloat *prob_hs)
{

	// read probability and define the host dynamic

        int tm = int((Tmax-Tnb[0])/(Tnb[1]-Tnb[0]));
	myfloat Ib[2][nB];

	for(int t = 0; t<2; t++ )
	{
		for(int nb = 0; nb<nB; nb++ )
		{
			Ib[t][nb] = I_blackC(Tnb[tm+t],narrBand[nb].wvc);
		}
	}

	for (int i=0; i<2; i++)
	{
		prob_hs[idx_p(i,0)] =  Ib[i][0]*(narrBand[0].wvr-narrBand[0].wvl)*100;
		for (int j=1; j<nB; j++)
		{
			prob_hs[idx_p(i,j)] = prob_hs[idx_p(i,j-1)] + Ib[i][j]*(narrBand[j].wvr-narrBand[j].wvl)*100;
		}
	}
	for (int i=0; i<2; i++)
	{
		for (int j=0; j<nB; j++)
		{
			prob_hs[idx_p(i,j)] = prob_hs[idx_p(i,j)]/prob_hs[idx_p(i,nB-1)];

		}
	}
}
