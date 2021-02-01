/*
 * definitions.h
 *
 *  Created on: May 8, 2017
 *      Author: simone
 */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_
#include "../param/NarrowBand.h"


struct beam {
        int ic, jc, kc;
        int i,  j,  k ;
        myfloat xp, yp ,zp;
        myfloat sx, sy, sz;
        myfloat Ti;
        myfloat tra;
};
typedef struct beam Beam;

struct Emission {
	myfloat west, east, north, south, top, bot;
};

struct EmissSpec {
	myfloat *west, *east;
	myfloat *north, *south;
	myfloat *top, *bot;
	void allocVar(int nb_)
	{
		west  = new myfloat[nb_];
		east  = new myfloat[nb_];
		north = new myfloat[nb_];
		south = new myfloat[nb_];
		top   = new myfloat[nb_];
		bot   = new myfloat[nb_];
	}
	void destroyVar()
	{
		delete[] west;
		delete[] east;
		delete[] north;
		delete[] south;
		delete[] top;
		delete[] bot;
	}
};

//Count structure
class Count
{
public:
	int nb_cnt[nB];
	int g_cnt[nQ][nB];
};

//Grid structure
class Gridn
{
public:
	myfloat *x,*xu,*xu0;
	myfloat *y,*yu;
	myfloat *z,*zu;
	int im,jm,km,im0;
	int sm;
	void mk_grid(int im_, int jm_, int km_, int im0_, int sm_)
	{

		jm = jm_;
		km = km_;
		im = im_;
		im0= im0_;
		sm= sm_;

		x  = new myfloat[im+2];
		xu = new myfloat[im+1];
		xu0 = new myfloat[im0+1];
		y  = new myfloat[jm+2];
		z  = new myfloat[km+2];
		yu = new myfloat[jm+1];
		zu = new myfloat[km+1];

		myfloat dy = Ly *1.0f / jm;
		myfloat dz = Lz *1.0f / km;

		for (int i = 0 ; i < jm+1 ; i++ )
			yu[i] = i * dy;
		for (int i = 0 ; i < km+1 ; i++ )
			zu[i] = i * dz;

		y[0] = -dy / 2;
		z[0] = -dz / 2;

		for (int i = 1 ; i < jm+1 ; i++ )
		{
			y[i] = (yu[i-1] + yu[i]) / 2;
		}
		for (int i = 1 ; i < km+1 ; i++ )
		{
			z[i] = (zu[i-1] + zu[i]) / 2;
		}
		y[jm+1] = y[jm]+dy;
		z[km+1] = z[km]+dz;


		xu0[0] = 0;
		for (int i=1; i < im0/2+1; i++)
		{
			myfloat r = (1.f*i)/im0;
			myfloat dx= 0.5 - fact_mesh*(r-0.5)*(r-0.5);
			xu0[i]=xu0[i-1]+dx;
		}
		myfloat xnorm=xu0[im0/2];
		for (int i=1; i < im0/2+1; i++)
		{
			xu0[i] = Lx/2.*xu0[i]/xnorm;
		}
		for (int i = im0; i>im0/2; i--)
		{
			xu0[i] = Lx - xu0[im0-i];
		}

		xu[0] = 0;
		for (int i=1; i < im+1; i++)
		{
			xu[i]=xu0[i*im0/im];
		}
		for( int i=1; i < im+1; i++)
			x[i] = 0.5*(xu[i]+xu[i-1]);

		x[0] = - x[1];
		x[im+1] = xu[im] + (xu[im] - x[im]);

	}
	void allocVar(int im_, int jm_, int km_, int sm_)
	{
		im = im_;
		jm = jm_;
		km = km_;
		sm = sm_;
		x  = new myfloat[im+2];
		y  = new myfloat[jm+2];
		z  = new myfloat[km+2];
		xu = new myfloat[im+1];
		yu = new myfloat[jm+1];
		zu = new myfloat[km+1];
	}
	void destroyVar()
	{
		delete[] x;
		delete[] y;
		delete[] z;
		delete[] xu;
		delete[] xu0;
		delete[] yu;
		delete[] zu;
	}
};
//Narrowband structure
class NarrowBand
{
public:
	myfloat kq[nT][nQ];
	myfloat wvc, wvl, wvr;
	myfloat kavg;
	int idx;
};

class Var_CPU
{
public:
	myfloat *x,*xu,*xu0;
	myfloat *y,*yu;
	myfloat *z,*zu;
	myfloat *T,*Tp,*C;
	int im,jm,km,im0;

	void mk_grid(int im_, int jm_, int km_, int im0_)
	{

		jm = jm_;
		km = km_;
		im = im_;
		im0= im0_;

		x  = new myfloat[im+2];
		xu = new myfloat[im+1];
		xu0 = new myfloat[im0+1];
		y  = new myfloat[jm+2];
		z  = new myfloat[km+2];
		yu = new myfloat[jm+1];
		zu = new myfloat[km+1];

		myfloat dy = Ly / jm;
		myfloat dz = Lz / km;

		for (int i = 0 ; i < jm+1 ; i++ )
			yu[i] = i * dy;
		for (int i = 0 ; i < km+1 ; i++ )
			zu[i] = i * dz;

		y[0] = -dy / 2;
		z[0] = -dz / 2;

		for (int i = 1 ; i < jm+1 ; i++ )
		{
			y[i] = (yu[i-1] + yu[i]) / 2;
		}
		for (int i = 1 ; i < km+1 ; i++ )
		{
			z[i] = (zu[i-1] + zu[i]) / 2;
		}
		y[jm+1] = y[jm]+dy;
		z[km+1] = z[km]+dz;


		xu0[0] = 0;
		for (int i=1; i < im0/2+1; i++)
		{
			myfloat r = (1.f*i)/im0;
			myfloat dx= 0.5 - 1.45*(r-0.5)*(r-0.5);
			xu0[i]=xu0[i-1]+dx;
		}
		myfloat xnorm=xu0[im0/2];
		for (int i=1; i < im0/2+1; i++)
		{
			xu0[i] = Lx/2.*xu0[i]/xnorm;
		}
		for (int i = im0; i>im0/2; i--)
		{
			xu0[i] = Lx - xu0[im0-i];
		}

		xu[0] = 0;
		for (int i=1; i < im+1; i++)
		{
			xu[i]=xu0[i*im0/im];
		}
		for( int i=1; i < im+1; i++)
			x[i] = 0.5*(xu[i]+xu[i-1]);

		x[0] = - x[1];
		x[im+1] = xu[im] + (xu[im] - x[im]);

		T   = new myfloat[(km+2)*(jm+2)*(im+2)];

	}
	void destroyVar()
	{
		delete[] x;
		delete[] y;
		delete[] z;
		delete[] xu;
		delete[] xu0;
		delete[] yu;
		delete[] zu;
		delete[] T;
	}
};

#endif /* DEFINITIONS_H_ */
