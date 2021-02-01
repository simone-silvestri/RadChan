
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <deque>

#define	 pi	 3.14159265358979323846
#define	 stefan	 5.670373e-08

#define nQuad 16
#define toll 1e-6
#define wvtoll 25
#define dwv 0.01
#define Deltawv1 50
#define Deltawv2 100
#define Deltawv3 100

#define n_pwrk 10000
#define nT 52
#define DT 10
#define Tinit 500

#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define pow3(a) \
   ({ __typeof__ (a) _a = (a); \
       _a*_a*_a; })

using namespace std;

class SolarBand
{
public:
  deque<double> iv, wv;
  double integ;
  double *prob;
  int size;
  void integral()
  {
	  integ=0;
	  prob = new double[iv.size()];
	  prob[0] =  (wv[1]/2)*iv[0];

	  for (int i=1; i<iv.size()-1; i++)
	  {
		  integ = integ+(wv[i+1]-wv[i-1])/2*(iv[i]) ;
		  prob[i] = integ;
	  }
	  integ = integ+(wv[iv.size()-1]-wv[iv.size()-2])*(iv[iv.size()-1]);
	  prob[iv.size()-1] = integ;
	  for(int i=0; i<iv.size()+1; i++)
	  {
		  prob[i] /= integ;
	  }
  }
};


class PartSpec2
{
public:
	int tot;
	char file[38];
	double *kv,*sv,*wv;
	void alloc_size()
	{
		kv = new double[tot];
		sv = new double[tot];
		wv = new double[tot];
		sprintf(file,"particle-data/carbon.txt");
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

double I_black( double T, double nu);
void readFillP(PartSpec2 *particles);
double lininterp(double y1, double y2, double x1, double x2, double x);

int main()
{
	char namefile[] = "solar-data/solarspec.txt";
	FILE *fp = fopen(namefile,"r");
	SolarBand solar_nm;
	SolarBand solar_cm;

	double iv,wv,dummy1;
	while (fscanf(fp,"%le %le %le %le \n",&wv,&dummy1,&dummy1,&iv) != EOF)
	{
    	solar_nm.iv.push_back(iv);  // convert from W/m2/nm in W/m2*cm
    	solar_nm.wv.push_back(wv);  // convert wavenumber from nm to in cm**(-1)
    	solar_cm.iv.push_front(iv*wv*wv/1e7);  // convert from W/m2/nm in W/m2*m
    	// taking into account the difference between a band in wavelengths and wavenumber
    	solar_cm.wv.push_front(1e7/wv);  // convert wavenumber from nm to in cm**(-1)
	}
	fclose(fp);

	solar_nm.integral();
	solar_cm.integral();
	solar_cm.size = solar_cm.wv.size();

	printf("the total solar radiation is: %f\n",solar_nm.integ);
	printf("the total solar radiation is: %f\n",solar_cm.integ);
	printf("the difference is is: %f\n",solar_cm.integ/solar_nm.integ);

	/****** CALCULATE THE PROBABILILTY DENSITY FUNCTIONS *******/

	solar_nm.prob[0] = 0;
	solar_cm.prob[0] = 0;

	/****** CALCULATE THE PARTICLE BANDS *******/

	PartSpec2 particles;
	particles.tot = solar_cm.size;
	particles.alloc_size();

	for(int i=0; i<solar_cm.size; i++)
		particles.wv[i] = solar_cm.wv[i];

	readFillP(&particles);

	/****** WRITE THE RESULTS *******/

	fp = fopen("tables_solar/prob.txt","w+");
	for (int i = 0; i<solar_cm.iv.size(); i++)
	{
		fprintf(fp,"%d\t%le\t%le\t%le\t%le\t%le\t%le\n",i,solar_cm.wv[i],solar_cm.iv[i],solar_cm.prob[i],solar_cm.integ,particles.kv[i],particles.sv[i]);
	}
	fclose(fp);

	fp = fopen("solar_blackbcm","w+");
	for (int i = 0; i<solar_cm.iv.size(); i++)
	{
		fprintf(fp,"%le\t%le\t%le\n",solar_cm.wv[i],solar_cm.iv[i],I_black(5782,solar_cm.wv[i]));
	}
	fclose(fp);


	return 0;
}

double I_black( double T, double nu)
{
    double C1 = 3.741771790075259e-16;
    double C2 = 0.014387741858429;
	return 1.0 / pi * C1 * pow3(nu*100) / ( expf((float)(C2*nu*100/T))-1);
}

void readFillP(PartSpec2 *particles)
{
	double wv,ev,sv,dummy1,dummy2;
	int dummy;
	int idx,h;
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
	for(int i = 0; i < particles->tot; i++)
	{
		// find index of partspec where wv lies and linear interpolation
		if(particles->wv[i] < partspec.wv[0])
		{
			particles->kv[i] = 0;
			particles->sv[i] = 0;
		}
		else if(particles->wv[i] > partspec.wv[partspec.wv.size()-1])
		{
			particles->kv[i] = partspec.kv[partspec.wv.size()-1];
			particles->sv[i] = partspec.sv[partspec.wv.size()-1];
		}
		else
		{
			for(h = 1; h < partspec.wv.size(); h++)
			{
				if((particles->wv[i] <= partspec.wv[h]) && (particles->wv[i] > partspec.wv[h-1])) {
					idx = h;
					break;
				}
			}
			particles->kv[i] = lininterp(partspec.kv[idx-1],partspec.kv[idx],partspec.wv[idx-1],partspec.wv[idx],particles->wv[i]);
			particles->sv[i] = lininterp(partspec.sv[idx-1],partspec.sv[idx],partspec.wv[idx-1],partspec.wv[idx],particles->wv[i]);
		}

	}

}

double lininterp(double y1, double y2, double x1, double x2, double x)
{
	return (y2 - y1) / (x2 - x1) * (x - x1) + y1;
}
