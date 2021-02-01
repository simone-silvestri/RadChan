/*
 * functions.h
 *
 *  Created on: Apr 26, 2017
 *      Author: simone
 */
#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_


//forward function definitions
void interp3D(Var_CPU *g);
void spline(myfloat x[], myfloat y[], int n, myfloat yp1, myfloat ypn, myfloat y2[]);
myfloat splint(myfloat xa[], myfloat ya[], myfloat y2a[], int n, myfloat x);
myfloat *vector12(long nl, long nh);
void free_vector(myfloat *v, long nl, long nh);

void readT(NarrowBand *narrBand, myfloat *Tnb, myfloat *kP , myfloat *prob_h, myfloat *probg_h, myfloat Tmax, myfloat *kappamax);
void readS(NarrowBand *narrBand, myfloat *Tnb, myfloat Tmax, myfloat *prob_hs);

#endif /* FUNCTIONS_H_ */
