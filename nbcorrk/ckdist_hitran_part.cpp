
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <deque>

#include "param.h"
#include "functions.h"

using namespace std;

/****** START OF THE MAIN PROGRAM ******/

int main()
{
	//remove stuff in the tables folder
	printf("\n");
	printf("removing tables\n");
	system("exec rm tables/*.txt");
	printf("\n");
	printf("removing parameters\n");
	system("exec rm errors/*");
	printf("\n");


	double T = 2500.;
	
	char namefile[]="FOLDER/2500K.txt";

	deque<PartSpecBand> bands;
	printf("namefile %s\n",namefile);
	
	printf("\n");
	double Entot = absWideBand(bands, namefile, T);

	printf("\n");
	FILE *fp = fopen("test.txt","wt");
	for (int i=0; i<bands.size(); i++)
	{
		bands[i].printBand(fp,i);
	}
	fclose(fp);

	double total = 0.0;
	for (int i=0; i<bands.size(); i++)
	{
		bands[i].alloc_ibv();
		for (int j=0; j<bands[i].kv.size(); j++)
			bands[i].ibv[j] = I_black(T, bands[i].wv[j] );

		bands[i].integral();
		total = total + bands[i].integ;
	}

	/******* STARTING TO REDUCE IN NARROW BANDS, BASED ON BAND NUMBER ********/


	// first we have to decide the number of narrow bands in every macro-band
	int totbands= 0;
	for (int i=0; i<bands.size(); i++)
	{
		if(i < 3)
		{
			bands[i].numbands =  (bands[i].wv[bands[i].wv.size()-1] - bands[i].wv[0])/Deltawv1 ;
			if(bands[i].numbands == 0)
				bands[i].numbands = 1;
			bands[i].sizebands = (bands[i].wv[bands[i].wv.size()-1] - bands[i].wv[0])/(double) (bands[i].numbands) ;
		}
		else if (i >= 3 && i < 6)
		{
			bands[i].numbands =  (bands[i].wv[bands[i].wv.size()-1] - bands[i].wv[0])/Deltawv2 ;
			if(bands[i].numbands == 0)
				bands[i].numbands = 1;
			bands[i].sizebands = (bands[i].wv[bands[i].wv.size()-1] - bands[i].wv[0])/(double) (bands[i].numbands) ;
		}
		else
		{
			bands[i].numbands =  (bands[i].wv[bands[i].wv.size()-1] - bands[i].wv[0])/Deltawv3 ;
			if(bands[i].numbands == 0)
				bands[i].numbands = 1;
			bands[i].sizebands = (bands[i].wv[bands[i].wv.size()-1] - bands[i].wv[0])/(double) (bands[i].numbands) ;
		}
	}
	for (int i=0; i<bands.size(); i++)
	{
		printf("Wide-band %-2d,\t Narrow-bands %-3d, \t delta wv %-2lf,\t Wide-band size %-2lf\n",i,bands[i].numbands,bands[i].sizebands,bands[i].wv[bands[i].wv.size()-1] - bands[i].wv[0]);
		totbands = totbands + bands[i].numbands;
	}
	printf("\n");
	printf("Total narrow bands : %d\n",totbands);
	printf("\n");

	/******* FILL NARROW BANDS WITH WAVENUMBERS AND  ABSORPTION COEFFICIENT ********/

	// allocate narrow bands
	for (int i=0; i<bands.size(); i++)
	{
		bands[i].alloc_narrbands();
		int h = 0;
		int comp = bands[i].kv.size()/bands[i].numbands;
		for (int j=0; j<bands[i].wv.size(); j++)
		{
			bands[i].narrband[h].kvn.push_back(bands[i].kv[j]);
			bands[i].narrband[h].wvn.push_back(bands[i].wv[j]);
			if ( j % comp == 0 && j != 0 && h != bands[i].numbands-1 )
				h=h+1;
		}
		for (int j=0; j<bands[i].numbands; j++)
		{
			bands[i].narrband[j].wv_left = bands[i].narrband[j].wvn[0];
			bands[i].narrband[j].wv_right = bands[i].narrband[j].wvn[bands[i].narrband[j].wvn.size()-1];
			bands[i].narrband[j].wvc = (bands[i].narrband[j].wv_left + bands[i].narrband[j].wv_right)/2;
			bands[i].narrband[j].ibv_left = I_black(T,bands[i].narrband[j].wv_left);
			bands[i].narrband[j].ibv_right = I_black(T,bands[i].narrband[j].wv_right);
			bands[i].narrband[j].ibvc = I_black(T,bands[i].narrband[j].wvc);
		}
	}

	fp = fopen("test2.txt","wt");
	for (int i=0; i<bands.size(); i++)
	{
		bands[i].print_narrbands(fp);
	}
	fclose(fp);

	/******** START WITH CORRELATED - K APPROACH ************/

	correlatedK(bands, totbands);

	/***************************************************************************************/
	/****** TEMPERATURE and PARTICLE ABSORPTION SPECTRA, READ, DISTRIBUTE AND WRITE ********/
	/***************************************************************************************/

	// allocate dimensions to the variables
	// write down probability table

	TempBand tempBand[nT];
//	TempBand nanoBand[nT];
#if Part==1
	PartSpec2 particles;
	PartSpec2 nanoPart;
	Phase phi;
#endif


	allocT(tempBand  , bands, totbands);
//	allocT(nanoBand  , bands, totbands);
#if Part==1
	allocP(&particles, bands, totbands);
//	allocN(&nanoPart , bands, totbands);
	allocH(&phi      , bands, totbands);
#endif

	// read the variables and fill the narrowBands

	readFillT(tempBand  , bands);
//	readFillT(nanoBand  , bands);
#if Part==1
	readFillP(&particles, bands);
//	readFillP(&nanoPart , bands);
	readFillH(&phi      , bands);
#endif

	//creating the nanofluid
//	for( int t=0; t<nT; t++)
//		for( int i=0; i<nanoPart.totbands; i++)
//			for( int j=0; j<nanoBand[t].narrband[i].size; j++ )
//			{
//				nanoBand[t].narrband[i].kvn[j] += nanoPart.narr_abs[i].kvn[j];
//			}

	// to test the interpolation of the spectrum un-comment the following line
	//particles.printSpec("interp-test");
	// to test the interpolation of the phase function un-comment the following line
//	phi.printPhi("interp-test");

	// correlated k for all the temperature files and particles

	correlatedKT(tempBand);
//	correlatedKT(nanoBand);
#if Part==1
	correlatedKP(&particles);
#endif 

	printf("\n");
	printf("Finished correlating k values\n");
	printf("\n");

	//calculating kP from the narrow bands to see if it is correct and probability
	//and calculating the probability of scattering angle

	calcProbT(tempBand);
//	calcProbT(nanoBand);
#if Part==1
	calcProbP(&particles);
	calcProbH(&phi);
#endif

	//uncomment the next line to test the angles probability
	//phi.printProb("angles-test");

	printf("\n");
	printf("Writing down tables\n");
	printf("\n");

	writeT(tempBand);
	printf("\n");
	printf("Finished writing gas tables\n");
	printf("\n");
//	writeN(nanoBand);

#if Part==1
	//if narrow band correlated k particle discretization is required uncomment the following line
	//writeP(&particles);
	//if narrow band box model particle discretization is required uncomment the following line
	writePbox(&particles);
	writeH(&phi);
#endif

	printf("\n");
	printf("Finished writing\n");
	printf("\n");

	return 0;

}


double I_black( double T, double nu)
{
	//myfloat k = 1.38065156e23;
	//myfloat h = 6.626070040e-34;
	//myfloat c0 = 299792458;
	//myfloat C1, C2;
	//C1 = 2 * pi * h * c0 * c0;
	//C2 = h * c0 / k;
	double C1 = 3.741771790075259e-16;
	double C2 = 0.014387741858429;

	return 1.0 / pi * C1 * pow((nu*100), 3) / ( expf((float)(C2*nu*100/T)) -1);

}


double I_deriv( double T, double nu)
{
	double C1 = 3.741771790075259e-16;
	double C2 = 0.014387741858429;

        C1 = C1* pow((nu*100), 3) / pi;
	C2 = C2* nu * 100;
        double C3 = expf((float)(C2/T));

	return (C1*C2*C3)/(T*T*(C3-1)*(C3-1));

}
