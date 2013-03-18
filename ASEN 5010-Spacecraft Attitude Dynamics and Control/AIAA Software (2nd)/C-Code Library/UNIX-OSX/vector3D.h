/*
 *  vector3D.h
 *  
 *
 *  Created by Hanspeter Schaub on Sat Mar 01 2003.
 *  Copyright (c) 2003 Hanspeter Schaub. All rights reserved.
 *
 */
#include <stdio.h>
#include <math.h>


#ifndef VECTOR3D
#define VECTOR3D

#ifdef __cplusplus
extern "C"  {
#endif



void 	set(double, double, double, double *);
void 	setMatrix(double, double, double, double, double, double, double, double, double, double mat[4][4]);
void 	set4(double, double, double, double , double *);
void  setZero(double *);
double 	norm(double *);
double 	norm4(double *);
void	cross(double *, double *, double *);
void	mult(double, double *, double *);
void  	equal(double *, double *);
double 	dot(double *, double *);
void	add(double *, double *, double *);
void	sub(double *, double *, double *);
void    printVector(const char *str, double *vec);
void    printVector4(const char *str, double *vec);
void    printMatrix(const char *str, double mat[4][4]);
int	    equalCheck(double *, double *, double d);
int     MequalCheck(double mat[4][4], double m2[4][4], double d);
int     VequalCheck(double *a, double *b, double d);
int     V4equalCheck(double *a, double *b, double d);

void	dotT(double *, double *, double mat[4][4]);
void	Mdot(double m[4][4], double *, double *);
void	MdotM(double m1[4][4], double m2[4][4], double ans[4][4]);
void	MdotMT(double m1[4][4], double m2[4][4], double ans[4][4]);
void	transpose(double m[4][4], double ans[4][4]);
void	Mmult(double, double m[4][4], double m2[4][4]);
void	Madd(double m1[4][4], double m2[4][4], double m3[4][4]);
void	Msub(double m1[4][4], double m2[4][4], double m3[4][4]); 
void	tilde(double *, double ans[4][4]);
void 	inverse(double m[4][4], double ans[4][4]);
double  detM(double a[4][4]);
double	trace(double m[4][4]);
void    eye(double m[4][4]);


double	arc_cosh(double);
double  arc_tanh(double);



/*
  common variables
*/
#ifndef M_PI
#define M_PI	    3.141592653589793
#endif

#ifndef D2R
#define D2R         M_PI/180.
#endif
#ifndef R2D
#define R2D         180./M_PI
#endif
#ifndef NAN
#define NAN         sqrt((float)-1)
#endif




#ifdef __cplusplus
}
#endif


#endif
