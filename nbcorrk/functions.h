/*
 * functions.h
 *
 *  Created on: May 4, 2017
 *      Author: simone
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

double absWideBand(deque<PartSpecBand> &bands, char *namefile, double T);
void allocT(TempBand *tempBand, deque<PartSpecBand> &bands, int totbands);
void allocP(PartSpec2 *partBand, deque<PartSpecBand> &bands, int totbands);
void allocN(PartSpec2 *particles, deque<PartSpecBand> &bands, int totbands);
void allocH(Phase *phi, deque<PartSpecBand> &bands, int totbands);
void readFillT(TempBand *tempBand, deque<PartSpecBand> &bands);
void readFillP(PartSpec2 *partBand, deque<PartSpecBand> &bands);
void readFillH(Phase *phi, deque<PartSpecBand> &bands);
void correlatedK(deque<PartSpecBand> &bands, int totbands);
void correlatedKT(TempBand *tempBand);
void correlatedKP(PartSpec2 *particles);
void calcProbT(TempBand *tempBand);
void calcProbP(PartSpec2 *particles);
void calcProbH(Phase *phi);
void writeT(TempBand *tempBand);
void writeP(PartSpec2 *particles);
void writeN(TempBand *tempBand);
void writePbox(PartSpec2 *particles);
void writeH(Phase *phi);
double I_black(double T, double nu);
double I_deriv(double T, double nu);

#endif /* FUNCTIONS_H_ */
