/*
 * classes.h
 *
 *  Created on: May 4, 2017
 *      Author: simone
 */

#ifndef CLASSES_H_
#define CLASSES_H_

#include <math.h>
#include <deque>

#define	 pi	 3.14159265358979323846
#define	 stefan	 5.670373e-08

#define nQuad 16 
#define toll 1e-6
#define wvtoll 25000
#define dwv 0.01
#define Deltawv1 10 
#define Deltawv2 10
#define Deltawv3 10

#define n_pwrk 20000
#define nT 2 //52 
#define nA 180
#define DT 10
#define Tinit 2500 
#define LBL   0 
#define H     1.0f
#define Part  0

#define FOLDER CO2-data

#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

using namespace std;

#if nQuad==100
	static const double gq[nQuad] = { 0.000140, 0.000739, 0.001816, 0.003370, 0.005400, 0.007904, 0.010879, 0.014322,
		0.018231, 0.022601, 0.027429, 0.032709, 0.038437, 0.044607, 0.051213, 0.058249, 0.065709, 0.073584, 0.081868,
		0.090552, 0.099629, 0.109089, 0.118924, 0.129123, 0.139678, 0.150578, 0.161813, 0.173372, 0.185243, 0.197417,
		0.209880, 0.222620, 0.235627, 0.248887, 0.262387, 0.276115, 0.290058, 0.304202, 0.318533, 0.333038, 0.347703,
		0.362514, 0.377456, 0.392516, 0.407679, 0.422930, 0.438255, 0.453640, 0.469068, 0.484527, 0.500000, 0.515473,
		0.530932, 0.546360, 0.561745, 0.577070, 0.592321, 0.607484, 0.622544, 0.637486, 0.652297, 0.666962, 0.681467,
		0.695798, 0.709942, 0.723885, 0.737613, 0.751113, 0.764373, 0.777380, 0.790120, 0.802583, 0.814757, 0.826628,
		0.838187, 0.849422, 0.860322, 0.870877, 0.881076, 0.890911, 0.900371, 0.909448, 0.918132, 0.926416, 0.934291,
		0.941751, 0.948787, 0.955393, 0.961563, 0.967291, 0.972571, 0.977399, 0.981769, 0.985678, 0.989121, 0.992096,
		0.994600, 0.996630, 0.998184, 0.999261};
	static const double wq[nQuad] = { 0.000360, 0.000838, 0.001316, 0.001792, 0.002267, 0.002740, 0.003210, 0.003677, 0.004140,
		0.004600, 0.005055, 0.005505, 0.005950, 0.006389, 0.006822, 0.007249, 0.007669, 0.008081, 0.008485, 0.008882, 0.009270,
		0.009649, 0.010019, 0.010379, 0.010729, 0.011069, 0.011399, 0.011717, 0.012024, 0.012320, 0.012604, 0.012876, 0.013135,
		0.013382, 0.013616, 0.013837, 0.014045, 0.014240, 0.014420, 0.014587, 0.014740, 0.014879, 0.015004, 0.015114, 0.015210,
		0.015291, 0.015357, 0.015409, 0.015446, 0.015468, 0.015476, 0.015468, 0.015446, 0.015409, 0.015357, 0.015291, 0.015210,
		0.015114, 0.015004, 0.014879, 0.014740, 0.014587, 0.014420, 0.014240, 0.014045, 0.013837, 0.013616, 0.013382, 0.013135,
		0.012876, 0.012604, 0.012320, 0.012024, 0.011717, 0.011399, 0.011069, 0.010729, 0.010379, 0.010019, 0.009649, 0.009270,
		0.008882, 0.008485, 0.008081, 0.007669, 0.007249, 0.006822, 0.006389, 0.005950, 0.005505, 0.005055, 0.004600, 0.004140,
		0.003677, 0.003210, 0.002740, 0.002267, 0.001792, 0.001316, 0.000838};
#elif nQuad==60
	static const double gq[nQuad] = { 0.000382, 0.002013, 0.004942, 0.009162, 0.014663, 0.021430, 0.029446,
		0.038689, 0.049135, 0.060758, 0.073527, 0.087409, 0.102367, 0.118362, 0.135353, 0.153295, 0.172142,
		0.191844, 0.212350, 0.233607, 0.255558, 0.278147, 0.301315, 0.325002, 0.349145, 0.373681, 0.398547,
		0.423678, 0.449008, 0.474471, 0.500000, 0.525529, 0.550992, 0.576322, 0.601453, 0.626319, 0.650855,
		0.674998, 0.698685, 0.721853, 0.744442, 0.766393, 0.787650, 0.808156, 0.827858, 0.846705, 0.864647,
		0.881638, 0.897633, 0.912591, 0.926473, 0.939242, 0.950865, 0.961311, 0.970554, 0.978570, 0.985337,
		0.990838, 0.995058, 0.997987};
	static const double wq[nQuad] = { 0.000981, 0.002280, 0.003576, 0.004863, 0.006137, 0.007395, 0.008633,
		0.009849, 0.011040, 0.012201, 0.013331, 0.014426, 0.015483, 0.016500, 0.017474, 0.018403, 0.019283,
		0.020113, 0.020890, 0.021613, 0.022280, 0.022889, 0.023437, 0.023925, 0.024350, 0.024712, 0.025009,
		0.025241, 0.025407, 0.025507, 0.025541, 0.025507, 0.025407, 0.025241, 0.025009, 0.024712, 0.024350,
		0.023925, 0.023437, 0.022889, 0.022280, 0.021613, 0.020890, 0.020113, 0.019283, 0.018403, 0.017474,
		0.016500, 0.015483, 0.014426, 0.013331, 0.012201, 0.011040, 0.009849, 0.008633, 0.007395, 0.006137,
		0.004863, 0.003576, 0.002280};
#elif nQuad==40
	static const double gq[nQuad] = { 0.000839, 0.004416, 0.010831, 0.020047, 0.032012, 0.046657, 0.063899,
		0.083639, 0.105764, 0.130148, 0.156649, 0.185118, 0.215390, 0.247292, 0.280641, 0.315247, 0.350912,
		0.387430, 0.424593, 0.462188, 0.500000, 0.537812, 0.575407, 0.612570, 0.649088, 0.684753, 0.719359,
		0.752708, 0.784610, 0.814882, 0.843351, 0.869852, 0.894236, 0.916361, 0.936101, 0.953343, 0.967988,
		0.979953, 0.989169, 0.995584};
	static const double wq[nQuad] = { 0.002153, 0.005000, 0.007822, 0.010601, 0.013318, 0.015959, 0.018509,
		0.020953, 0.023276, 0.025467, 0.027511, 0.029398, 0.031117, 0.032657, 0.034010, 0.035169, 0.036126,
		0.036876, 0.037415, 0.037739, 0.037848, 0.037739, 0.037415, 0.036876, 0.036126, 0.035169, 0.034010,
		0.032657, 0.031117, 0.029398, 0.027511, 0.025467, 0.023276, 0.020953, 0.018509, 0.015959, 0.013318,
		0.010601, 0.007822, 0.005000};
#elif nQuad==17
	static const double gq[nQuad] = {0.0000000, 0.15540585, 0.45000000, 0.74459415, 0.90000000, 0.91554058,
			0.94500000, 0.97445942, 0.99000000, 0.99155406, 0.99450000, 0.99744594, 0.99900000, 0.99917267,
			0.99950000, 0.99982733, 1.00000000};
	static const double wq[nQuad] = {0.045, 0.245, 0.32, 0.245, 0.0495, 0.0245, 0.032, 0.0245, 0.00495, 0.00245,
			0.0032, 0.00245, 0.0005, 0.0002722, 0.00035556, 0.000272222, 0.00005};
#elif nQuad==16
	static const double gq[nQuad] = { 0.04830766,0.14447196,0.23928736,0.33186860,0.42135127,0.50689990,
		0.58771575,0.66304426,0.73218211,0.79448379,0.84936761,0.89632115,0.93490607,0.96476225,0.98561151,
		0.99726386};
	static const double wq[nQuad] = { 0.09654009,0.09563872,0.09384440,0.09117387,0.08765209,0.08331193,
		0.07819389,0.07234580,0.06582222,0.05868409,0.05099806,0.04283589,0.03427387,0.02539207,0.01627439,
		0.00701862};
#elif nQuad==10
	static const double gq[nQuad] = { 0.010886, 0.056469, 0.134924, 0.240452, 0.365228, 0.500000, 0.634772, 0.759548, 0.865076, 0.943531};
	static const double wq[nQuad] = { 0.027834, 0.062790, 0.093145, 0.116597, 0.131402, 0.136463, 0.131402, 0.116597, 0.093145, 0.062790};
#elif nQuad==7
	static const double gq[nQuad] ={ 0.00000 ,0.15541 ,0.45000 ,0.74459 ,0.90000 ,0.93551 ,0.98449};
	static const double wq[nQuad] = {0.04500,0.24500,0.32000,0.24500,0.05611,0.05125,0.03764};
#endif

class NarrowBand
{
public:
   deque<double> wvn,kvn;
   double wv_left, wv_right, wvc;
   double ibv_left, ibv_right, ibvc;
   double kmin,kmax;
   double ff[n_pwrk+1];
   double gg[n_pwrk+1];
   double kk[n_pwrk+1];
   double kq[nQuad];
};

class NarrowBand2
{
public:
	double *wvn,*kvn,*prob_lbl,*cuprob_lbl;
	int size;
	double wv_left, wv_right, wvc;
    double ibv_left, ibv_right, ibvc;
    double kmin,kmax;
    double prob;
    double prob_new2[nQuad];
    double k_avg;
    double cuprob;
    double cuprob_new2[nQuad];
    double ff[n_pwrk+1];
    double gg[n_pwrk+1];
    double kk[n_pwrk+1];
    double kq[nQuad];
	char file[38];
    void allocVar(int i)
    {
	  wvn = new double[i];
	  kvn = new double[i];
	  prob_lbl= new double[i];
	  cuprob_lbl= new double[i];
	}

};

class NarrowBandP
{
public:
	double *wvn,*kvn,*prob_lbl,*cuprob_lbl;
	int size;
	double wv_left, wv_right, wvc;
    double kmin,kmax,k_avg;
    double prob[nT];
    double prob_new2[nT][nQuad];
    double cuprob[nT];
    double cuprob_new2[nT][nQuad];
    double ff[n_pwrk+1];
    double gg[n_pwrk+1];
    double kk[n_pwrk+1];
    double kq[nQuad];
	char file[38];
    void allocVar(int i)
    {
	  wvn = new double[i];
	  kvn = new double[i];
	  prob_lbl= new double[i];
	  cuprob_lbl= new double[i];
	}
};

class TempBand
{
public:
	NarrowBand2 *narrband;
	double T;
    double kP;
    double kR;
    double kM;
    double kPr;
	int totbands;
	char file[38];
	char file_err[38];
	char file_par[38];
	void allocNarrBands(int i)
	{
	     narrband = new NarrowBand2[i];
	}
	void print_narrbands(FILE *fp)
	{
		for (int i=0; i<totbands; i++)
		{
			for (int j=0; j<narrband[i].size; j++)
		    fprintf(fp,"%d\t%d\t%le\t%le\n",i,j,narrband[i].wvn[j],narrband[i].kvn[j]);
			i++;
		}
	}
	void print_tempbands_quad(FILE *fp)
	{
		for (int i=0; i<totbands; i++)
		{
			for (int g=0; g<nQuad; g++)
			{
				fprintf(fp,"%le %le\n",(narrband[i].wv_left+(narrband[i].wv_right-narrband[i].wv_left)*gq[g]),narrband[i].kq[g]);
			}
		}
	}
	void print_tempbands_totnb(FILE *fp)
	{
		for (int i=0; i<totbands; i++)
		{
			for (int g=0; g<n_pwrk+1; g++)
			{
				fprintf(fp,"%le %le\n",(narrband[i].wv_left+(narrband[i].wv_right-narrband[i].wv_left)*narrband[i].gg[g]),narrband[i].kk[g]);
			}
		}
	}
};

class PartSpecBand
{
public:
  deque<double> kv, wv;
  double *ibv;
  double integ;
  int numbands;
  double sizebands;

  NarrowBand *narrband;

  void printBand(FILE *fp, int j)
  {
     for (int i=0; i<kv.size(); i++)
        fprintf(fp,"%le\t%le\t%d\n",wv[i],kv[i],j);
  }
  void integral()
  {
  integ=0;
  for (int i=1; i<kv.size(); i++)
   integ = integ+(wv[i]-wv[i-1])*(kv[i]+kv[i-1])*(ibv[i]+ibv[i-1])/4;
  }
  void alloc_ibv()
  {
     ibv = new double[kv.size()];
  }
  void alloc_narrbands()
  {
     narrband = new NarrowBand[numbands];
  }
  void print_narrbands(FILE *fp)
  {
	 for (int i=0; i<numbands; i++)
	 {
		 for (int j=0; j<narrband[i].wvn.size(); j++)
	      fprintf(fp,"%d\t%d\t%le\t%le\n",i,j,narrband[i].wvn[j],narrband[i].kvn[j]);
		 i++;
	 }
  }
};

class PartSpec
{
public:
	deque<double> kv, wv, sv;
	void printSpec(const char *filename)
	{
		FILE *fpf = fopen(filename,"w+");
		for(int j = 0; j < kv.size(); j++)
		{
			fprintf(fpf,"%lf %lf %lf\n",wv[j],kv[j],sv[j]);
		}
		fclose(fpf);
	}
};

class PartSpec2
{
public:
	NarrowBandP *narr_abs;
	NarrowBandP *narr_sca;
	int totbands;
	double T[nT], kP[nT];
	char file[38];
	char file_err[38];
	char file_par[38];
	void allocNarrBands(int i)
	{
	     narr_abs = new NarrowBandP[i];
	     narr_sca = new NarrowBandP[i];
	}
	void printSpec(const char *filename)
	{
		FILE *fpf = fopen(filename,"w+");
		for(int i = 0; i < totbands; i++)
			for(int j = 0; j < narr_abs[i].size; j++)
			{
				fprintf(fpf,"%lf %lf %lf\n",narr_abs[i].wvn[j],narr_abs[i].kvn[j],narr_sca[i].kvn[j]);
			}
		fclose(fpf);
	}
	void calcPlanckabs()
	{
		T[0] = Tinit;
		for (int i = 1; i < nT; i++)
			T[i] = T[i-1] + DT;

		for (int i = 0; i < nT; i++)
		{
			kP[i] = 0;
			for (int j=0; j < totbands; j++)
			{
				for (int h=1; h<n_pwrk; h++)
				{
					kP[i] += narr_abs[j].kk[h]*I_part(T[i], narr_abs[j].wvc)*(narr_abs[j].ff[h])*
							(narr_abs[j].wv_right - narr_abs[j].wv_left)*100*pi/(stefan * T[i]*T[i]*T[i]*T[i]);
				}
			}
		}
	}
	double I_part( double T, double nu)
	{
		double C1 = 3.741771790075259e-16;
		double C2 = 0.014387741858429;
		return 1.0 / pi * C1 * pow((nu*100), 3) / ( expf((float)(C2*nu*100/T)) -1);
	}
};

class Phaset
{
public:
	deque<double> wv, phiv[nA];
	char file[50];
	double tv[nA];
	void printPhi(const char *filename)
	{
		FILE *fpf = fopen(filename,"w+");
		for(int i = 0; i < wv.size(); i++)
			for(int j = 0; j < nA; j++)
			{
				tv[j] = j*pi/(nA*1.0);
				fprintf(fpf,"%lf %lf %lf\n",wv[i],tv[j],phiv[j][i]);
			}
		fclose(fpf);
	}
};

class Phase
{
public:
	double *p[nA], t[nA], *wvc;
	int totbands;
	char file[50];
	char writefile[50];
	double *prob[nA], *cuprob[nA];
	void allocVar(size_t i)
	{
		wvc = new double[i];
		for(int h = 0; h < nA; h++)
		{
			p[h] = new double[i];
			prob[h] = new double[i];
			cuprob[h] = new double[i];
		}
	}
	void printPhi(const char *filename)
	{
		FILE *fpf = fopen(filename,"w+");
		for(int i = 0; i < totbands; i++)
			for(int j = 0; j < nA; j++)
			{
				fprintf(fpf,"%lf %lf %lf\n",wvc[i],t[j],p[j][i]);
			}
		fclose(fpf);
	}
	void printProb(const char *filename)
	{
		FILE *fpf = fopen(filename,"w+");
		for(int i = 0; i < totbands; i++)
			for(int j = 0; j < nA; j++)
			{
				fprintf(fpf,"%lf %lf %lf %lf\n",wvc[i],t[j],cuprob[j][i],prob[j][i]);
			}
		fclose(fpf);
	}
};

#endif /* CLASSES_H_ */
