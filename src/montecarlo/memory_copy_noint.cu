
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
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand_kernel.h>
#include "../param/param.h"
#include "definitions.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include "memory.h"

extern myfloat probg[nT][nB][nQ];

void black_copy(EmissSpec *Ibw_d, NarrowBand *narrBand)
{

	EmissSpec *Ibw_h;
	Ibw_h = (EmissSpec*)malloc(sizeof(EmissSpec));
	Ibw_h->allocVar(nB);

    for(int nb = 0; nb<nB; nb++ )
    {
#if srt ==1
    	Ibw_h->west[nb]  = I_blackC(Tww,narrBand[narrBand[nb].idx].wvc);
    	Ibw_h->east[nb]  = I_blackC(Twe,narrBand[narrBand[nb].idx].wvc);
    	Ibw_h->north[nb] = I_blackC(Twn,narrBand[narrBand[nb].idx].wvc);
    	Ibw_h->south[nb] = I_blackC(Tws,narrBand[narrBand[nb].idx].wvc);
    	Ibw_h->top[nb]   = I_blackC(Twt,narrBand[narrBand[nb].idx].wvc);
    	Ibw_h->bot[nb]   = I_blackC(Twb,narrBand[narrBand[nb].idx].wvc);
# else
    	Ibw_h->west[nb]  = I_blackC(Tww,narrBand[nb].wvc);
    	Ibw_h->east[nb]  = I_blackC(Twe,narrBand[nb].wvc);
    	Ibw_h->north[nb] = I_blackC(Twn,narrBand[nb].wvc);
    	Ibw_h->south[nb] = I_blackC(Tws,narrBand[nb].wvc);
    	Ibw_h->top[nb]   = I_blackC(Twt,narrBand[nb].wvc);
    	Ibw_h->bot[nb]   = I_blackC(Twb,narrBand[nb].wvc);
#endif
    }

	myfloat *west,*east;
	myfloat *north,*south;
	myfloat *top,*bot;
	EmissSpec *temp_Ib;

	temp_Ib = (EmissSpec*)malloc(nB*sizeof(EmissSpec));

	cudaMalloc((void**)&west , nB * sizeof(myfloat));
	cudaMalloc((void**)&east , nB * sizeof(myfloat));
	cudaMalloc((void**)&north, nB * sizeof(myfloat));
	cudaMalloc((void**)&south, nB * sizeof(myfloat));
	cudaMalloc((void**)&top  , nB * sizeof(myfloat));
	cudaMalloc((void**)&bot  , nB * sizeof(myfloat));

	cudaMemcpy(west , Ibw_h->west , nB * sizeof(myfloat), cudaMemcpyHostToDevice);
	cudaMemcpy(east , Ibw_h->east , nB * sizeof(myfloat), cudaMemcpyHostToDevice);
	cudaMemcpy(north, Ibw_h->north, nB * sizeof(myfloat), cudaMemcpyHostToDevice);
	cudaMemcpy(south, Ibw_h->south, nB * sizeof(myfloat), cudaMemcpyHostToDevice);
	cudaMemcpy(top  , Ibw_h->top  , nB * sizeof(myfloat), cudaMemcpyHostToDevice);
	cudaMemcpy(bot  , Ibw_h->bot  , nB * sizeof(myfloat), cudaMemcpyHostToDevice);

	temp_Ib->west  = west;
	temp_Ib->east  = east;
	temp_Ib->north = north;
	temp_Ib->south = south;
	temp_Ib->top   = top;
	temp_Ib->bot   = bot;

	cudaMemcpy(Ibw_d, temp_Ib, nB * sizeof(EmissSpec), cudaMemcpyHostToDevice);
	cudaCheckErrors("cudaMalloc Ibw error");

	Ibw_h->destroyVar();
	free(Ibw_h);
        free(temp_Ib);
}


void grid_copy(Gridn *gridCPU, Gridn *gridGPU)
{
    myfloat *xGPU[grid_num],*xuGPU[grid_num];
    myfloat *yGPU[grid_num],*yuGPU[grid_num];
    myfloat *zGPU[grid_num],*zuGPU[grid_num];
    Gridn *temp_G;
    temp_G = (Gridn*)malloc(grid_num*sizeof(Gridn));

    for (int grd=0; grd<grid_num; grd++)
    {
    	cudaMalloc((void**)&xGPU[grd] , (gridCPU[grd].im+2) * sizeof(myfloat));
    	cudaMalloc((void**)&xuGPU[grd], (gridCPU[grd].im+1) * sizeof(myfloat));
    	cudaMalloc((void**)&yGPU[grd] , (gridCPU[grd].jm+2) * sizeof(myfloat));
    	cudaMalloc((void**)&yuGPU[grd], (gridCPU[grd].jm+1) * sizeof(myfloat));
    	cudaMalloc((void**)&zGPU[grd] , (gridCPU[grd].km+2) * sizeof(myfloat));
    	cudaMalloc((void**)&zuGPU[grd], (gridCPU[grd].km+1) * sizeof(myfloat));

    	cudaMemcpy(xGPU[grd] , gridCPU[grd].x , (gridCPU[grd].im+2) * sizeof(myfloat), cudaMemcpyHostToDevice);
    	cudaMemcpy(xuGPU[grd], gridCPU[grd].xu, (gridCPU[grd].im+1) * sizeof(myfloat), cudaMemcpyHostToDevice);
    	cudaMemcpy(yGPU[grd] , gridCPU[grd].y , (gridCPU[grd].jm+2) * sizeof(myfloat), cudaMemcpyHostToDevice);
    	cudaMemcpy(yuGPU[grd], gridCPU[grd].yu, (gridCPU[grd].jm+1) * sizeof(myfloat), cudaMemcpyHostToDevice);
    	cudaMemcpy(zGPU[grd] , gridCPU[grd].z , (gridCPU[grd].km+2) * sizeof(myfloat), cudaMemcpyHostToDevice);
    	cudaMemcpy(zuGPU[grd], gridCPU[grd].zu, (gridCPU[grd].km+1) * sizeof(myfloat), cudaMemcpyHostToDevice);

    	temp_G[grd].x  = xGPU[grd];
    	temp_G[grd].xu = xuGPU[grd];
    	temp_G[grd].y  = yGPU[grd];
    	temp_G[grd].yu = yuGPU[grd];
    	temp_G[grd].z  = zGPU[grd];
    	temp_G[grd].zu = zuGPU[grd];
    	temp_G[grd].im = gridCPU[grd].im;
    	temp_G[grd].jm = gridCPU[grd].jm;
    	temp_G[grd].km = gridCPU[grd].km;
    	temp_G[grd].sm = gridCPU[grd].sm;

    }
	cudaMemcpy(gridGPU, temp_G, grid_num * sizeof(Gridn), cudaMemcpyHostToDevice);
	cudaCheckErrors("cudaMalloc Grid error");
        free(temp_G);
}

void temp_fluid_copy(cudaTextureObject_t *tex_tempf_d, Var_CPU *varCPU)
{
     cudaArray *T_volumeArray;
     static cudaTextureObject_t tex_temp;
     
     //First grid temperature
     const cudaExtent extentT = make_cudaExtent(varCPU[0].im+2, varCPU[0].jm+2, varCPU[0].km+2);
     cudaChannelFormatDesc channelDescT = cudaCreateChannelDesc<myfloat>();
     cudaMalloc3DArray(&T_volumeArray, &channelDescT, extentT);
     cudaCheckErrors("cudaMalloc3D error");
     
     // Copying host memory to 3D cuda_Array
     cudaMemcpy3DParms copyParamsT = {0};
     copyParamsT.srcPtr   = make_cudaPitchedPtr((void*)varCPU[0].T, extentT.width*sizeof(myfloat), extentT.width, extentT.height);
     copyParamsT.dstArray = T_volumeArray;
     copyParamsT.extent   = extentT;
     copyParamsT.kind     = cudaMemcpyHostToDevice;
     
     cudaMemcpy3D(&copyParamsT);
     cudaCheckErrors("cudaMemcpy3D fail");
     
     // Binding 3D cuda_Array to texture array
     cudaResourceDesc    texRes;
     memset(&texRes, 0, sizeof(cudaResourceDesc));
     texRes.resType = cudaResourceTypeArray;
     texRes.res.array.array  = T_volumeArray;
     cudaTextureDesc     texDescr;
     memset(&texDescr, 0, sizeof(cudaTextureDesc));
     texDescr.normalizedCoords = false;
     texDescr.filterMode = cudaFilterModeLinear;
     texDescr.addressMode[0] = cudaAddressModeClamp;   // clamp
     texDescr.addressMode[1] = cudaAddressModeClamp;
     texDescr.addressMode[2] = cudaAddressModeClamp;
     texDescr.readMode = cudaReadModeElementType;
     cudaCreateTextureObject(&tex_temp, &texRes, &texDescr, NULL);
     cudaCheckErrors("Bind fail");
     cudaMemcpy(tex_tempf_d , &tex_temp , sizeof(cudaTextureObject_t) ,cudaMemcpyHostToDevice);
}

void narrowband_copy(NarrowBand *narrBand, myfloat *wvc_d, int *idx_d, cudaTextureObject_t *tex_d, cudaTextureObject_t *tex_prob_d,
					myfloat *Tnb_d, myfloat *Tnb)
{
    myfloat *wvc_h;
    wvc_h = (myfloat*)malloc(nB*sizeof(myfloat));
    int *idx_h;
    idx_h = (int*)malloc(nB*sizeof(int));
    for(int nb = 0; nb<nB; nb++)
    	idx_h[nb] = narrBand[nb].idx;

    for(int nb = 0; nb<nB; nb++ )
    {
#if srt == 1
    	wvc_h[nb] = narrBand[narrBand[nb].idx].wvc;
#else
    	wvc_h[nb] = narrBand[nb].wvc;
#endif
    }

    cudaMemcpy(idx_d   , idx_h   , nB *      sizeof(int)     ,cudaMemcpyHostToDevice);
    cudaMemcpy(wvc_d   , wvc_h   , nB *      sizeof(myfloat) ,cudaMemcpyHostToDevice);
    cudaMemcpy(Tnb_d   , Tnb     , nT *      sizeof(myfloat) ,cudaMemcpyHostToDevice);

    free(idx_h);
    free(wvc_h);

    myfloat kq[nQ][nT][nB];

    for(int g = 0; g<nQ; g++)
    {
    	for(int t = 0; t<nT; t++)
		{
    		for(int nb = 0; nb<nB; nb++)
    		{
#if srt == 1
    			kq[g][t][nb] = narrBand[narrBand[nb].idx].kq[t][g];
#else
    			kq[g][t][nb] = narrBand[nb].kq[t][g];
#endif
    		}
		}
    }

    cudaArray *d_volumeArray = 0;
    static cudaTextureObject_t tex1_d;

    cudaArray *PG_volumeArray = 0;
	static cudaTextureObject_t tex1_prob_d;

	const cudaExtent extent = make_cudaExtent(nB, nT, nQ);
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<myfloat>();
	cudaMalloc3DArray(&d_volumeArray, &channelDesc, extent);
	cudaCheckErrors("cudaMalloc3D error");

	// Copying host memory to 3D cuda_Array
	cudaMemcpy3DParms copyParams = {0};
	copyParams.srcPtr   = make_cudaPitchedPtr((void*)kq, extent.width*sizeof(myfloat), extent.width, extent.height);
	copyParams.dstArray = d_volumeArray;
	copyParams.extent   = extent;
	copyParams.kind     = cudaMemcpyHostToDevice;

	cudaMemcpy3D(&copyParams);
	cudaCheckErrors("cudaMemcpy3D fail");

	// Binding 3D cuda_Array to texture array
	cudaResourceDesc    texRes;
	memset(&texRes, 0, sizeof(cudaResourceDesc));
	texRes.resType = cudaResourceTypeArray;
	texRes.res.array.array  = d_volumeArray;
	cudaTextureDesc     texDescr;
	memset(&texDescr, 0, sizeof(cudaTextureDesc));
	texDescr.normalizedCoords = false;
	texDescr.filterMode = cudaFilterModeLinear;
	texDescr.addressMode[0] = cudaAddressModeClamp;   // clamp
	texDescr.addressMode[1] = cudaAddressModeClamp;
	texDescr.addressMode[2] = cudaAddressModeClamp;
	texDescr.readMode = cudaReadModeElementType;
    cudaCreateTextureObject(&tex1_d, &texRes, &texDescr, NULL);
    cudaCheckErrors("Bind fail");

    cudaMemcpy(tex_d , &tex1_d , sizeof(cudaTextureObject_t) ,cudaMemcpyHostToDevice);


	const cudaExtent extentP = make_cudaExtent(nQ, nB, nT);
	cudaChannelFormatDesc channelDescP = cudaCreateChannelDesc<myfloat>();
	cudaMalloc3DArray(&PG_volumeArray, &channelDescP, extentP);
	cudaCheckErrors("cudaMalloc3D error");

	// Copying host memory to 3D cuda_Array
	cudaMemcpy3DParms copyParamsP = {0};
	copyParamsP.srcPtr   = make_cudaPitchedPtr((void*)probg, extentP.width*sizeof(myfloat), extentP.width, extentP.height);
	copyParamsP.dstArray = PG_volumeArray;
	copyParamsP.extent   = extentP;
	copyParamsP.kind     = cudaMemcpyHostToDevice;

	cudaMemcpy3D(&copyParamsP);
	cudaCheckErrors("cudaMemcpy3D fail");

	// Binding 3D cuda_Array to texture array
	cudaResourceDesc    texResP;
	memset(&texResP, 0, sizeof(cudaResourceDesc));
	texResP.resType = cudaResourceTypeArray;
	texResP.res.array.array  = PG_volumeArray;
	cudaTextureDesc     texDescrP;
	memset(&texDescrP, 0, sizeof(cudaTextureDesc));
	texDescrP.normalizedCoords = false;
	texDescrP.filterMode = cudaFilterModeLinear;
	texDescrP.addressMode[0] = cudaAddressModeClamp;   // clamp
	texDescrP.addressMode[1] = cudaAddressModeClamp;
	texDescrP.addressMode[2] = cudaAddressModeClamp;
	texDescrP.readMode = cudaReadModeElementType;
    cudaCreateTextureObject(&tex1_prob_d, &texResP, &texDescrP, NULL);
    cudaCheckErrors("Bind fail");

    cudaMemcpy(tex_prob_d , &tex1_prob_d , sizeof(cudaTextureObject_t) ,cudaMemcpyHostToDevice);

}

myfloat I_blackC( myfloat T, myfloat nu)
{
	// way of calculating C1 and C2 (nu is in cm^(-1) while the output is in W/m^2

    myfloat C1 = 3.741771790075259e-16;
    myfloat C2 = 0.014387741858429;

	return 1.0 / pi * C1 * pow3(nu*100) / (expf(C2*nu*100/T)-1);

}
