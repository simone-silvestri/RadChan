
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
#include <string.h> //for strcat

#include "param.h"
#include "functions.h"

void writeT(TempBand *tempBand)
{

	//write down plank-mean tables

	FILE *fp = fopen("../properties/planck-mean.txt","wt");
	for (int i=0; i<nT; i++)
		fprintf(fp,"%lf %le %le %le\n",tempBand[i].T,tempBand[i].kPr,tempBand[i].kR,tempBand[i].kM);
	fclose(fp);

	//write down probability table CK
	fp = fopen("../properties/prob.txt","wt");
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			fprintf(fp,"%lf \t %-5d \t%30.20le \n",tempBand[i].T,j,tempBand[i].narrband[j].cuprob);
		}
	}
	fclose(fp);

	//write down probability table CK new - 2
	fp = fopen("../properties/prob_new2.txt","wt");
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			fprintf(fp,"%lf \t %-5d \t%30.20le \t",tempBand[i].T,j,tempBand[i].narrband[j].cuprob);
			for (int g=0; g<nQuad; g++)
			{
				fprintf(fp,"%30.20le\t",tempBand[i].narrband[j].cuprob_new2[g]);
			}
			fprintf(fp,"\n");
		}
	}
	fclose(fp);
#if LBL == 1
	//write down probability table LBL
	fp = fopen("../properties/wave_lbl.txt","wt");
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			for (int g=0; g<tempBand[i].narrband[j].size; g++)
			{
				fprintf(fp,"%le \n",tempBand[i].narrband[j].wvn[g]);
			}
		}
	}
	fclose(fp);
	fp = fopen("../properties/prob_lbl.txt","wt");
	for (int i=0; i<nT; i++)
	{
		int n=0;
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			for (int g=0; g<tempBand[i].narrband[j].size; g++)
			{
				fprintf(fp,"%lf \t %-5d \t %30.20le \t %30.20le\n",tempBand[i].T,n,tempBand[i].narrband[j].cuprob_lbl[g],tempBand[i].narrband[j].kvn[g]);
				n=n+1;
			}
		}
		printf("Total elements in temperature %f is %d\n",tempBand[i].T,n);
	}
	fclose(fp);
#endif

	//write down tables of nb-T-kq;
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<tempBand[i].totbands; j++)
			sprintf(tempBand[i].narrband[j].file,"../properties/NarrBand%d.txt",j);
	}
	for (int j=0; j<tempBand[0].totbands; j++)
	{
		fp = fopen(tempBand[0].narrband[j].file,"wt");
		fprintf(fp,"%le\t%le\t%le \n",tempBand[0].narrband[j].wv_left,tempBand[0].narrband[j].wv_right,tempBand[0].narrband[j].wvc);
		for (int i=0; i<nT; i++)
		{
			fprintf(fp," %lf\t",tempBand[i].T);
			fprintf(fp," %le\t",tempBand[i].narrband[j].cuprob);
			fprintf(fp," %le\t",tempBand[i].narrband[j].k_avg);
			for (int h=0; h<nQuad; h++)
			{
				fprintf(fp," %le\t",tempBand[i].narrband[j].kq[h]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
	}

}

void writeN(TempBand *tempBand)
{

	//write down plank-mean tables

	FILE *fp = fopen("../properties/nanofluid/planck-mean.txt","wt");
	for (int i=0; i<nT; i++)
		fprintf(fp,"%lf %le %le\n",tempBand[i].T,tempBand[i].kPr,tempBand[i].kP);
	fclose(fp);

	//write down probability table CK
	fp = fopen("../properties/nanofluid/prob.txt","wt");
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			fprintf(fp,"%lf \t %-5d \t%30.20le \n",tempBand[i].T,j,tempBand[i].narrband[j].cuprob);
		}
	}
	fclose(fp);

	//write down probability table CK new - 2
	fp = fopen("../properties/nanofluid/prob_new2.txt","wt");
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			fprintf(fp,"%lf \t %-5d \t%30.20le \t",tempBand[i].T,j,tempBand[i].narrband[j].cuprob);
			for (int g=0; g<nQuad; g++)
			{
				fprintf(fp,"%30.20le\t",tempBand[i].narrband[j].cuprob_new2[g]);
			}
			fprintf(fp,"\n");
		}
	}
	fclose(fp);

#if LBL == 1
	//write down probability table LBL
	fp = fopen("../properties/nanofluid/wave_lbl.txt","wt");
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			for (int g=0; g<tempBand[i].narrband[j].size; g++)
			{
				fprintf(fp,"%le \n",tempBand[i].narrband[j].wvn[g]);
			}
		}
	}
	fclose(fp);
	fp = fopen("../properties/nanofluid/prob_lbl.txt","wt");
	for (int i=0; i<nT; i++)
	{
		int n=0;
		for (int j=0; j<tempBand[i].totbands; j++)
		{
			for (int g=0; g<tempBand[i].narrband[j].size; g++)
			{
				fprintf(fp,"%lf \t %-5d \t %30.20le \t %30.20le\n",tempBand[i].T,n,tempBand[i].narrband[j].cuprob_lbl[g],tempBand[i].narrband[j].kvn[g]);
				n=n+1;
			}
		}
		printf("Total elements in temperature %f is %d\n",tempBand[i].T,n);
	}
	fclose(fp);
#endif 

	//write down tables of nb-T-kq;
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<tempBand[i].totbands; j++)
			sprintf(tempBand[i].narrband[j].file,"../properties/nanofluid/NarrBand%d.txt",j);
	}
	for (int j=0; j<tempBand[0].totbands; j++)
	{
		fp = fopen(tempBand[0].narrband[j].file,"wt");
		fprintf(fp,"%le\t%le\t%le \n",tempBand[0].narrband[j].wv_left,tempBand[0].narrband[j].wv_right,tempBand[0].narrband[j].wvc);
		for (int i=0; i<nT; i++)
		{
			fprintf(fp," %lf\t",tempBand[i].T);
			fprintf(fp," %le\t",tempBand[i].narrband[j].cuprob);
			fprintf(fp," %le\t",tempBand[i].narrband[j].k_avg);
			for (int h=0; h<nQuad; h++)
			{
				fprintf(fp," %le\t",tempBand[i].narrband[j].kq[h]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
	}

}

void writeP(PartSpec2 *particles)
{
	//write down plank-mean tables
	FILE *fp = fopen("../properties/particles/planck-mean.txt","wt");
	for (int i=0; i<nT; i++)
		fprintf(fp,"%lf %le\n",particles->T[i],particles->kP[i]);
	fclose(fp);

	//write down probability table CK
	fp = fopen("../properties/particles/prob.txt","wt");
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<particles->totbands; j++)
		{
			fprintf(fp,"%lf \t %-5d \t%30.20le \n",particles->T[i],j,particles->narr_abs[j].cuprob[i]);
		}
	}
	fclose(fp);

	//write down probability table CK new - 2
	fp = fopen("../properties/particles/prob_new2.txt","wt");
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<particles->totbands; j++)
		{
			fprintf(fp,"%lf \t %-5d \t%30.20le \t",particles->T[i],j,particles->narr_abs[j].cuprob[i]);
			for (int g=0; g<nQuad; g++)
			{
				fprintf(fp,"%30.20le\t",particles->narr_abs[j].cuprob_new2[i][g]);
			}
			fprintf(fp,"\n");
		}
	}
	fclose(fp);

	//write down tables of nb-kq;
	for (int j=0; j < particles->totbands; j++)
		sprintf(particles->narr_abs[j].file,"../properties/particles/NarrBand%d.txt",j);

	for (int j=0; j < particles->totbands; j++)
	{
		fp = fopen(particles->narr_abs[j].file,"wt");
		fprintf(fp,"%le\t%le\t%le \n",particles->narr_abs[j].wv_left,particles->narr_abs[j].wv_right,particles->narr_abs[j].wvc);
		for (int h=0; h<nQuad; h++)
		{
			fprintf(fp," %le\t",particles->narr_abs[j].kq[h]);
		}
		fprintf(fp,"\n");
		for (int h=0; h<nQuad; h++)
		{
			fprintf(fp," %le\t",particles->narr_sca[j].kq[h]);
		}
		fclose(fp);
	}

}

void writePbox(PartSpec2 *particles)
{
	//write down plank-mean tables
	FILE *fp = fopen("../properties/particles/planck-mean.txt","wt");
	for (int i=0; i<nT; i++)
		fprintf(fp,"%lf %le\n",particles->T[i],particles->kP[i]);
	fclose(fp);

	//write down probability table CK
	fp = fopen("../properties/particles/prob.txt","wt");
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<particles->totbands; j++)
		{
			fprintf(fp,"%lf \t %-5d \t%30.20le \n",particles->T[i],j,particles->narr_abs[j].cuprob[i]);
		}
	}
	fclose(fp);

#if LBL == 1
	//write down probability table LBL
	fp = fopen("../properties/particles/wave_lbl.txt","wt");
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<particles->totbands; j++)
		{
			for (int g=0; g<particles->narr_abs[j].size; g++)
			{
				fprintf(fp,"%le \n",particles->narr_abs[j].wvn[g]);
			}
		}
	}
	fclose(fp);
	fp = fopen("../properties/particles/prob_lbl.txt","wt");
	for (int i=0; i<nT; i++)
	{
		int n=0;
		for (int j=0; j<particles->totbands; j++)
		{
			for (int g=0; g<particles->narr_abs[j].size; g++)
			{
				fprintf(fp,"%lf \t %-5d \t %30.20le \t %30.20le\n",particles->T[i],n,particles->narr_abs[j].cuprob_lbl[g],particles->narr_abs[j].kvn[g]);
				n=n+1;
			}
		}
		printf("Total elements in temperature %f is %d\n",particles->T[i],n);
	}
	fclose(fp);
#endif
	//write down tables of nb-kq;
	for (int j=0; j < particles->totbands; j++)
		sprintf(particles->narr_abs[j].file,"../properties/particles/NarrBand%d.txt",j);

	for (int j=0; j < particles->totbands; j++)
	{
		fp = fopen(particles->narr_abs[j].file,"wt");
		fprintf(fp,"%le\t%le\t%le\t%le\t%le \n",particles->narr_abs[j].wv_left,particles->narr_abs[j].wv_right,particles->narr_abs[j].wvc,
									  particles->narr_abs[j].k_avg,  particles->narr_abs[j].k_avg);
		fclose(fp);
	}

}


void writeH(Phase *phi)
{
	sprintf(phi->writefile,"../properties/particles/prob_scatt.txt");
	FILE *fp = fopen(phi->writefile,"w+");
	for(int i = 0; i < phi->totbands; i++)
		for(int j = 0; j < nA; j++)
		{
			fprintf(fp,"%le\t%le\t%le\n",phi->t[j],phi->wvc[i],phi->cuprob[j][i]);
		}
	fclose(fp);
}
