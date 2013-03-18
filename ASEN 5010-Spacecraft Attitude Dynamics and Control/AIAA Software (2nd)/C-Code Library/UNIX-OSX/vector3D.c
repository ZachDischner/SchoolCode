/*
 *  vector3D.c
 *  
 *
 *  Created by Hanspeter Schaub on Sat Mar 01 2003.
 *  Copyright (c) 2003 Hanspeter Schaub. All rights reserved.
 *
 *
 *	Provides a library for doing 3D matrix algebra.
 *
 */

#include "vector3D.h"



void    set(double x, double y, double z, double *v) 
{
    v[1] = x;
    v[2] = y;
    v[3] = z;
}

void 	setMatrix(double m11, double m12, double m13, 
                double m21, double m22, double m23, 
                double m31, double m32, double m33, double mat[4][4])
{
  mat[1][1] = m11;  mat[1][2] = m12;  mat[1][3] = m13;
  mat[2][1] = m21;  mat[2][2] = m22;  mat[2][3] = m23;
  mat[3][1] = m31;  mat[3][2] = m32;  mat[3][3] = m33;
}

void    set4(double x, double y, double z, double zz, double *v) 
{
    v[1] = x;
    v[2] = y;
    v[3] = z;
    v[4] = zz;
}

void    setZero(double *v)
{
  v[1] = 0.;
  v[2] = 0.;
  v[3] = 0.;
}

double	norm(double *x)
{
	double ans;
	
	ans = sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
	
	return ans;
}

double	norm4(double *x)
{
	double ans;
	
	ans = sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]);
	
	return ans;
}

void	cross(double *x, double *y, double *h)
{
	h[1] = x[2]*y[3] - x[3]*y[2];
	h[2] = x[3]*y[1] - x[1]*y[3];
	h[3] = x[1]*y[2] - x[2]*y[1];
}


void	mult(double g, double *x, double *ans)
{
	ans[1] = g*x[1];
	ans[2] = g*x[2];
	ans[3] = g*x[3];
}

void	equal(double *y, double *x)
{
    x[1] = y[1];
    x[2] = y[2];
    x[3] = y[3];
	
	return;
}

double	dot(double *x, double *y)
{
	double ans;
	
	ans = x[1]*y[1] + x[2]*y[2] + x[3]*y[3];
	
	return ans;
}


void	add(double *x, double *y, double *ans)
{
	ans[1] = x[1] + y[1];
	ans[2] = x[2] + y[2];
	ans[3] = x[3] + y[3];
}

void	sub(double *x, double *y, double *ans)
{
	ans[1] = x[1] - y[1];
	ans[2] = x[2] - y[2];
	ans[3] = x[3] - y[3];
}









void    printVector(const char *str, double *vec)
{
    printf("%s (%g, %g, %g)\n", str, vec[1], vec[2], vec[3]);
     
    return;
}

void    printVector4(const char *str, double *vec)
{
    printf("%s (%g, %g, %g, %g)\n", str, vec[1], vec[2], vec[3], vec[4]);
     
    return;
}

void    printMatrix(const char *str, double mat[4][4])
{
    printf("%s:\n\t\t%12.8g\t%12.8g\t%12.8g\n", str, mat[1][1], mat[1][2], mat[1][3]);
    printf("\t\t%12.8g\t%12.8g\t%12.8g\n",         mat[2][1], mat[2][2], mat[2][3]);
    printf("\t\t%12.8g\t%12.8g\t%12.8g\n",         mat[3][1], mat[3][2], mat[3][3]);
     
    return;
}

int	    equalCheck(double *v1, double *v2, double d)
{
  int c = 1, i;
  
  for (i=1;i<=3;i++) {
	if (fabs(v1[i] - v2[i]) > pow(10.,-d))
	  c = 0;
  }
  
  return c;
}


int     MequalCheck(double m1[4][4], double m2[4][4], double d)
{
  int c = 1, i,j;
  
  for (i=1;i<=3;i++) {
	for (j=1;j<=3;j++) {
	  if (fabs(m1[i][j] - m2[i][j]) > pow(10.,-d)) {
		c = 0;
	  }
	}
  }
  
  return c;
}


int     VequalCheck(double *a, double *b, double d)
{
  int c = 1, i;
  
  for (i=1;i<=3;i++) {
	if (fabs(a[i] - b[i]) > pow(10.,-d))
	  c = 0;
  }
  
  return c;
}

int     V4equalCheck(double *a, double *b, double d)
{
  int c = 1, i;
  
  for (i=1;i<=4;i++) {
	if (fabs(a[i] - b[i]) > pow(10.,-d))
	  c = 0;
  }
  
  return c;
}






void	dotT(double *v1, double *v2, double ans[4][4]) {
  mult(v1[1], v2, ans[1]);
  mult(v1[2], v2, ans[2]);
  mult(v1[3], v2, ans[3]);
  
  return;
}





void	Mdot(double m[4][4], double *a, double *b)
{
  b[1] = dot(m[1],a);
  b[2] = dot(m[2],a);
  b[3] = dot(m[3],a);
}

void    transpose(double m1[4][4], double ans[4][4])
{
  ans[1][1]=m1[1][1];  ans[1][2]=m1[2][1];  ans[1][3]=m1[3][1];
  ans[2][1]=m1[1][2];  ans[2][2]=m1[2][2];  ans[2][3]=m1[3][2];
  ans[3][1]=m1[1][3];  ans[3][2]=m1[2][3];  ans[3][3]=m1[3][3];
  
  return;
}


void	MdotM(double m1[4][4], double m0[4][4], double ans[4][4])
{
  double m2[4][4];
  transpose(m0,m2);
  
  ans[1][1] = dot(m1[1],m2[1]);
  ans[1][2] = dot(m1[1],m2[2]);
  ans[1][3] = dot(m1[1],m2[3]);
  ans[2][1] = dot(m1[2],m2[1]);
  ans[2][2] = dot(m1[2],m2[2]);
  ans[2][3] = dot(m1[2],m2[3]);
  ans[3][1] = dot(m1[3],m2[1]);
  ans[3][2] = dot(m1[3],m2[2]);
  ans[3][3] = dot(m1[3],m2[3]);
  
  return;
}




void	MdotMT(double m1[4][4], double m0[4][4], double ans[4][4])
{
  
  ans[1][1] = dot(m1[1],m0[1]);
  ans[1][2] = dot(m1[1],m0[2]);
  ans[1][3] = dot(m1[1],m0[3]);
  ans[2][1] = dot(m1[2],m0[1]);
  ans[2][2] = dot(m1[2],m0[2]);
  ans[2][3] = dot(m1[2],m0[3]);
  ans[3][1] = dot(m1[3],m0[1]);
  ans[3][2] = dot(m1[3],m0[2]);
  ans[3][3] = dot(m1[3],m0[3]);
  
  return;
}






void	Mmult(double a, double m[4][4], double m2[4][4])
{
  mult(a, m[1], m2[1]);
  mult(a, m[2], m2[2]);
  mult(a, m[3], m2[3]);
  
  return;
}


void	Madd(double m1[4][4], double m2[4][4], double m3[4][4]) {
  add(m1[1], m2[1], m3[1]);
  add(m1[2], m2[2], m3[2]);
  add(m1[3], m2[3], m3[3]);
  
  return;
}


void	Msub(double m1[4][4], double m2[4][4], double m3[4][4]) {
  sub(m1[1], m2[1], m3[1]);
  sub(m1[2], m2[2], m3[2]);
  sub(m1[3], m2[3], m3[3]);
  
  return;
}



void	tilde(double *a, double ans[4][4])
{
  ans[1][1] = ans[2][2] = ans[3][3] = 0.;
  ans[1][2] = -a[3];
  ans[2][1] =  a[3];
  ans[1][3] =  a[2];
  ans[3][1] = -a[2];
  ans[2][3] = -a[1];
  ans[3][2] =  a[1];
  
  return;
} 


double detM(double m[4][4])
{
  return -m[1][3]*m[2][2]*m[3][1] + m[1][2]*m[2][3]*m[3][1] + m[1][3]*m[2][1]*m[3][2]
		 -m[1][1]*m[2][3]*m[3][2] - m[1][2]*m[2][1]*m[3][3] + m[1][1]*m[2][2]*m[3][3];
}



/*
	finds the inverse of a (3 x 3) sized matrix
*/
void 	inverse(double m[4][4], double ans[4][4])
{
	double	det;
	
	det = detM(m);
    
	if (fabs(det) > 1e-12) {
	  ans[1][1] = (-m[2][3]*m[3][2] + m[2][2]*m[3][3])/det;
	  ans[1][2] = (m[1][3]*m[3][2]-m[1][2]*m[3][3])/det;
	  ans[1][3] = (-m[1][3]*m[2][2]+m[1][2]*m[2][3])/det;
	  ans[2][1] = (m[2][3]*m[3][1]-m[2][1]*m[3][3])/det;
	  ans[2][2] = (-m[1][3]*m[3][1]+m[1][1]*m[3][3])/det;
	  ans[2][3] = (m[1][3]*m[2][1]-m[1][1]*m[2][3])/det;
	  ans[3][1] = (-m[2][2]*m[3][1]+m[2][1]*m[3][2])/det;
	  ans[3][2] = (m[1][2]*m[3][1]-m[1][1]*m[3][2])/det;
	  ans[3][3] = (-m[1][2]*m[2][1]+m[1][1]*m[2][2])/det;
	} else {
	  fprintf(stderr, "ERROR: singular 3x3 matrix inverse\n");
	  ans[1][1]=ans[1][2]=ans[1][3]=ans[2][1]=ans[2][2]=ans[2][3]=
		ans[3][1]=ans[3][2]=ans[3][3]= NAN;
	}
	
	
	return;
}


double	trace(double m[4][4])
{
  return m[1][1] + m[2][2] + m[3][3];
}


void	eye(double m[4][4])
{
  m[1][1] = 1.0;	m[1][2] = 0.0;	  m[1][3] = 0.0;
  m[2][1] = 0.0;	m[2][2] = 1.0;	  m[2][3] = 0.0;
  m[3][1] = 0.0;	m[3][2] = 0.0;	  m[3][3] = 1.0;
  
  return;
}


double	arc_cosh(double z)
{
  return log(z+sqrt(z+1.)*sqrt(z-1.));
}

double arc_tanh(double z)
{
  return 0.5*(log(1.+z) - log(1.-z));
} 
