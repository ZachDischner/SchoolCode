/*
 *  orbitalMotion.h
 *  OrbitalMotion
 *
 *  Created by Hanspeter Schaub on 6/19/05.
 *  
 *  This package provides various orbital 
 *  mechanics subroutines using in astrodynamics calculations.
 *
 */

#include <stdio.h>
#include <math.h>
#include "astroConstants.h"
#include "vector3D.h"

#ifndef ORBITALMOTION
#define ORBITALMOTION

#ifdef __cplusplus
extern "C"  {
#endif


typedef struct classicElem {
	double a;
	double e;
	double i;
	double Omega;
	double omega;
	double anom;
} classicElements;


double	E2f(double E, double e);
double	E2M(double E, double e);
double	f2E(double f, double e);
double	f2H(double f, double e);
double	H2f(double H, double e);
double	H2N(double H, double e);
double	M2E(double M, double e);
double	N2H(double N, double e);
void	elem2rv(double mu, classicElements *elements, double *rVec, double *vVec);
void    rv2elem(double mu, double *rVec, double *vVec, classicElements *elements);

double AtmosphericDensity( double alt );
void   AtmosphericDrag(double Cd, double A, double m, double *rvec, double *vvec, double *advec);
void   JPerturb(double *rvec, int num, double *ajtot);


#define N_DEBYE_PARAMETERS 37


#ifdef __cplusplus
}
#endif


#endif
