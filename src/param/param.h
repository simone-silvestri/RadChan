/*******************************************************************//**
*     NARROWBAND MONTE CARLO OERM GPU by Simone Silvestri
*     instructions in doc.h
**********************************************************************/

//#include "doc.h"

#define  num_gpu    4 
#define  p_col      48       
#define  imax       406  
#define  jmax       576  
#define  kmax       2304  
#define  Lx         2.f
#define  Ly         6.283185307179586f
#define  Lz         37.69911184307751f   //62.83185307179586f 
#define  fact_mesh  1.6
//variance block calculator and photon per block
#define  nVar       1  
#define  photmax    6000 
#define  chanHeight 0.5


//number of grid used (maximum 4, triggers preprocessor #if statements)
#define grid_num 9 

// REMEMBER! If INTERP='y' then subsequent meshes should be divisible between each other.
const int maxi[] = {imax, imax  , imax  , imax   , imax   , imax/29, imax/29, imax/29 , imax/203};
const int maxj[] = {jmax, jmax/3, jmax/6, jmax/12, jmax/24, jmax/48, jmax/96, jmax/96 , jmax/192};
const int maxk[] = {kmax, kmax/2, kmax/4, kmax/8 , kmax/16, kmax/32, kmax/64, kmax/128, kmax/256};
//maximum number of steps in the grids
#if grid_num==1
const int maxs[] = {100000};
#elif grid_num==2
const int maxs[] = {5, 100000};
#elif grid_num==3
const int maxs[] = {5, 5, 100000};
#elif grid_num==4
const int maxs[] = {5, 5, 5, 100000};
#elif grid_num==5
const int maxs[] = {5, 5, 5, 5, 100000};
#elif grid_num==6
const int maxs[] = {5, 5, 5, 5, 5, 100000};
#elif grid_num==7
const int maxs[] = {5, 5, 5, 5, 5, 5, 100000};
#elif grid_num==8
const int maxs[] = {5, 5, 5, 5, 5, 5, 5, 100000};
#elif grid_num==9
const int maxs[] = {5, 5, 5, 5, 5, 5, 5, 5, 100000};
#endif

//GPU details
#define num_streams     2 
#define blockdim        128	
#define nblocks         4096	

//sorting or not sorting the narrow bands (1 -> yes, other int -> no)
#define srt 1
//grey radiation (1 -> yes, other int -> no)
#define grey 0 
//calculate walls (1 -> yes, other int -> no)
#define calcwalls 0 
//random sequence (1 -> sobol, other int -> pseudo-random)
#define random 1 
// 1 -> double precision, 0 -> single precision
#define floatpoint 0

//define the boundary values
#define	 epsw	 1.f
#define	 epse	 1.f
#define	 epss	 1.f
#define	 epsn	 1.f
#define	 epsb	 1.f
#define	 epst	 1.f

#define	 Tww	 0 // not important if adiabw==1
#define	 Twe	 0 // not important if adiabe==1
#define	 Tws	 0
#define	 Twn	 0
#define	 Twb	 0 
#define	 Twt	 0 

#define	 adiabw	 1 
#define	 adiabe	 1  

#define	 stefan	 5.670373e-08
#define	 toll	 1.0e-3
#define	 tollS	 1.0e-3
#define	 pi 	 3.14159265358979323846

// 1 -> periodic walls
// 2 -> black walls
// 3 -> specular reflection 
// 4 -> black walls with zero gradient temperature
// 5 -> lost ray (temperature equal to the emitter) 
#define bdw 2 
#define bde 2 
#define bds 1    
#define bdn 1 
#define bdb 4 
#define bdt 4 

#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)

#define cudaCheckErrors2(msg,rank) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d). We are in rank %d\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__, rank); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)

#define idx_p(t,b) \
		({ t*nB + b; }) //for values going from 0 to max-1, 2-D
#define idx_I(t,b) \
		({ t + b*nT; }) //for values going from 0 to max-1, 2-D
#define idx_pg(t,b,g) \
		({ g + nQ*( b + t*nB ) ; }) //for values going from 0 to max-1, 3-D
#define idx_T(i,j,k,im,jm) \
		({ i + (im+2)*( j + k*(jm+2) ) ; }) //for values going from 0 to max-1, 3-D
#define idx_S(i,k) \
		({ (k-1) + (i-1) * kmax ;} ) //for values going from 1 to max 2-D
#define idx_D(i,j,k) \
        ({ (k-1)+kmax*((j-1)+(i-1)*jmax); }) //for values going from 1 to max 3-D
#define idx_F(i,j,k) \
		({ i + (imax+2)*( j + k*(jmax+2) ) ; }) //for values going from 0 to max-1, 3-D

#define pow4(a) \
        ({ a*a*a*a ; })

#define pow3(a) \
        ({ a*a*a ; })

#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a >= _b ? _a : _b; })
#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a <= _b ? _a : _b; })

// modifying accordingly the read statements

#if floatpoint == 1
#define myfloat double
#define reade "%le"
#define readf "%lf"
#elif floatpoint == 0
#define myfloat float
#define reade "%e"
#define readf "%f"
#endif

#if adiabw==1
#define boundini 0
#else
#define boundini 1 
#endif
#if adiabe==1
#define boundend 2 
#else
#define boundend 1 
#endif

#define myfloatF float 
