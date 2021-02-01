/*
 * device_functions.h
 *
 *  Created on: May 8, 2017
 *      Author: simone
 */

#ifndef DEVICE_FUNCTIONS_H_
#define DEVICE_FUNCTIONS_H_

__global__ void kernel_fluid(curandDirectionVectors32_t *rngVectors, Gridn *my_grid, myfloat *wvc_d, myfloat *Tnb_t, myfloat *var,
					   myfloat *device, int n, int ns, int stream, myfloat kappamax, myfloat Tmax, int *idx_nb,
					   cudaTextureObject_t *tex_Tf, cudaTextureObject_t *tex, cudaTextureObject_t *tex_prob,
					   int gpu, int zstart, EmissSpec *Ibw);
__global__ void kernel_surface(Gridn *my_grid, myfloat *wvc_d, myfloat *Tnb_t, int n, int ns, int stream, myfloat Tmax, int *idx_nb,
		cudaTextureObject_t *tex_Tf, cudaTextureObject_t *tex, cudaTextureObject_t *tex_prob,
		int gpu, int zstart);
__device__ __forceinline__  void emiss_angSobol(Beam *ray, curandStateSobol32_t *state, curandStateSobol32_t *state2);
__device__ __forceinline__  void emiss_ang(Beam *ray, curandState_t *state);
__device__ __forceinline__  void emiss_angS(Beam *ray, curandState_t *state, int bound);
__device__ __forceinline__ void wave_find(int *nb, int *g, int *countnb, int *countg, int *idx_nb, Count *cnt);
__device__ __forceinline__ void find_band(myfloat *Tnb, myfloat Tmax, int *g, int *nb, curandState_t *state, int tm, cudaTextureObject_t *tex_prob);
__device__ __forceinline__ void find_bandS(myfloat *Tnb, myfloat Tmax, int *g, int *nb, curandState_t *state, int tm);
__device__ __forceinline__ void march_ray(Beam *ray, int nb, int g, int *flag, myfloat *De_OERMc,
						myfloat Ibmax, Gridn grid, myfloat ratioI, myfloat invTi4, cudaTextureObject_t tex_Tf,
						Emission emiss, myfloat wvc, const myfloat Tnb[], cudaTextureObject_t *tex, myfloat tollerance);
__device__ __forceinline__  void kernel_find(Count *count, myfloat *Tnb, myfloat Tmax, curandState_t *state, cudaTextureObject_t *tex_prob);
__device__ __forceinline__  void kernel_findS(Count *count, myfloat *Tnb, myfloat Tmax, curandState_t *state);
__global__ void kernel_results (myfloat *var, myfloat *device, myfloat *varS, myfloat *devS, int n, int ns, int stream);
//__global__ void kernel_resultsS(int ns, int stream);
#endif /* DEVICE_FUNCTIONS_H_ */
