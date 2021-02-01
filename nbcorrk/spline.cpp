
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

#include "nrutil.h"


class Temp
{
public:
	myfloat *T;
	int im,jm,km;
	void allocateT(int im_, int jm_, int km_)
	{
		im = im_;
		jm = jm_;
		km = km_;
		T  = new myfloat[(km+2)*(jm+2)*(im+2)];
	}
};

void interp3D(Var_CPU *g)
{
	Temp  tx[grid_num];
	Temp txy[grid_num];

	for(int grd = 1; grd < grid_num; grd++)
	{
		  tx[grd].allocateT(g[grd].im, jmax, kmax);
		 txy[grd].allocateT(g[grd].im, g[grd].jm, kmax);
	}

	//start the grid loop
	for(int j = 1; j < jmax+1; j++)
	{
		for(int k = 1; k < kmax+1; k++)
		{
			myfloat temp[imax+2],temp2[imax+2];
			for(int i = 1; i < imax+1; i++)
			{
				temp[i] = g[0].T[idx_T(i,j,k,imax,jmax)];
			}
			spline(g[0].x, temp, imax, 1.e30, 1.0e30, temp2 );
			for(int grd = 1; grd < grid_num; grd++)
			{
				for(int i = 1; i < tx[grd].im+1; i++)
				{
					tx[grd].T[idx_T(i,j,k,tx[grd].im,tx[grd].jm)] = splint(g[0].x,temp,temp2,imax,g[grd].x[i]);
				}
			}
		}
	}

	for(int grd = 1; grd < grid_num; grd++)
	{
		for(int i = 1; i < tx[grd].im+1; i++)
		{
			for(int k = 1; k < kmax+1; k++)
			{
				myfloat temp[jmax+2],temp2[jmax+2];
				for(int j = 1; j < jmax+1; j++)
				{
					temp[j] = tx[grd].T[idx_T(i,j,k,tx[grd].im,tx[grd].jm)];
				}
				spline(g[0].y, temp, jmax, 1e30, 1e30, temp2);
				for(int j = 1; j < txy[grd].jm+1; j++)
				{
					txy[grd].T[idx_T(i,j,k,txy[grd].im,txy[grd].jm)] = splint(g[0].y,temp,temp2,jmax,g[grd].y[j]);
				}
			}
		}

	}

	for(int grd = 1; grd < grid_num; grd++)
	{
		for(int i = 1; i < txy[grd].im+1; i++)
		{
			for(int j = 1; j < txy[grd].jm+1; j++)
			{
				myfloat temp[kmax+2],temp2[kmax+2];
				for(int k = 1; k < kmax+1; k++)
				{
					temp[k] = txy[grd].T[idx_T(i,j,k,txy[grd].im,txy[grd].jm)];
				}
				spline(g[0].z, temp, kmax, 1e30, 1e30, temp2);
				for(int k = 1; k < g[grd].km+1; k++)
				{
					g[grd].T[idx_T(i,j,k,g[grd].im,g[grd].jm)] = splint(g[0].z,temp,temp2,kmax,g[grd].z[k]);
				}
			}
		}

	}
	return;
}


void spline(myfloat x[], myfloat y[], int n, myfloat yp1, myfloat ypn, myfloat y2[])
{
	int i,j,k;
	myfloat p,qn,sig,un,*u;
	u=vector12(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_vector(u,1,n-1);
}

myfloat splint(myfloat xa[], myfloat ya[], myfloat y2a[], int n, myfloat x)
{
	int klo,khi,k;
	myfloat h,b,a,y;
	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0){
		fprintf(stderr,"Bad xa input to routine splint");
		fprintf(stderr,"...now exiting to system...\n");
		exit(1);
	}
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
	return y;
}

myfloat *vector12(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	myfloat *v;

	v=(myfloat *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(myfloat)));
	if (!v){
		fprintf(stderr,"allocation failure in vector()");
		fprintf(stderr,"...now exiting to system...\n");
		exit(1);
	}
	return v-nl+NR_END;
}
void free_vector(myfloat *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
