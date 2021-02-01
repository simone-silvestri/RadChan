/*
 * memory.h
 *
 *  Created on: May 9, 2017
 *      Author: simone
 */

#ifndef MEMORY_H_
#define MEMORY_H_

void black_copy(EmissSpec *Ibw_d, NarrowBand *narrBand);
void grid_copy(Gridn *gridCPU, Gridn *gridGPU);
void temp_fluid_copy(cudaTextureObject_t *tex_tempf_d, Var_CPU *varCPU);
void narrowband_copy(NarrowBand *narrBand, myfloat *wcv_d, int *idx_d,
		cudaTextureObject_t *tex_d, cudaTextureObject_t *tex_prob_d,
		myfloat *Tnb_d, myfloat *Tnb);
myfloat I_blackC( myfloat T, myfloat nu);

#endif /* MEMORY_H_ */
