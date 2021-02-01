/*******************************************************************//**
*
* NARROWBAND MONTE CARLO OERM GPU by Simone Silvestri
*           
*******************************************************************/

/** Using a narrow-band correlated-k discretization of the absorption spectra */
#if grey == 0             

#define nT 135            /*!< number of temperature bins in the narrow-band tables */      
#define nQ 16             /*!< number of quadrature points used for the correlated-k integration */
#define nB 999            /*!< number of narrow bands */ 

/** Gaussian quadrature points (gq) and weights (wq)  */
#if nQ==17
	static const myfloat gq[nQ] = {0.0000000, 0.15540585, 0.45000000, 0.74459415, 0.90000000, 0.91554058,
		0.94500000, 0.97445942, 0.99000000, 0.99155406, 0.99450000, 0.99744594, 0.99900000, 0.99917267,
		0.99950000, 0.99982733, 1.00000000};
	static const myfloat wq[nQ] = {0.045, 0.245, 0.32, 0.245, 0.0495, 0.0245, 0.032, 0.0245, 0.00495, 0.00245,
		0.0032, 0.00245, 0.0005, 0.0002722, 0.00035556, 0.000272222, 0.00005};
#elif nQ==16
	static const myfloat gq[nQ] = { 0.04830766,0.14447196,0.23928736,0.33186860,0.42135127,0.50689990,
		0.58771575,0.66304426,0.73218211,0.79448379,0.84936761,0.89632115,0.93490607,0.96476225,0.98561151,
		0.99726386};
	static const myfloat wq[nQ] = { 0.09654009,0.09563872,0.09384440,0.09117387,0.08765209,0.08331193,
		0.07819389,0.07234580,0.06582222,0.05868409,0.05099806,0.04283589,0.03427387,0.02539207,0.01627439,
		0.00701862};
#elif nQ==7
	static const myfloat gq[nQ] ={ 0.00000 ,0.15541 ,0.45000 ,0.74459 ,0.90000 ,0.93551 ,0.98449};
	static const myfloat wq[nQ] = {0.04500,0.24500,0.32000,0.24500,0.05611,0.05125,0.03764};
#endif        

/** Using a grey gas formulation for the participating media */
#else         
#define nT 65 
#define nQ 1 
#define nB 1000
#define abscoeff 1.f
#define Tmins 490.f
#define wvmax 5000.f
#undef srt
#define srt 0
#endif
