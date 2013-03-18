/*
 *  RigidBodyKinematics.c
 *  
 *
 *  Created by Hanspeter Schaub on Thu Nov 10 2005.
 *  Copyright (c) 2005 Hanspeter Schaub. All rights reserved.
 *
 *
 *	Provides a library for doing 3D matrix algebra.
 *
 */


#include "RigidBodyKinematics.h"




/*
% addEP(B1,B2)
%
%	Q = addEP(B1,B2) provides the Euler parameter vector
%	which corresponds to performing to successive
%	rotations B1 and B2.
%
*/
void addEP(double *b1, double *b2, double *q)
{
	
	q[1] = b2[1]*b1[1]-b2[2]*b1[2]-b2[3]*b1[3]-b2[4]*b1[4];
	q[2] = b2[2]*b1[1]+b2[1]*b1[2]+b2[4]*b1[3]-b2[3]*b1[4];
	q[3] = b2[3]*b1[1]-b2[4]*b1[2]+b2[1]*b1[3]+b2[2]*b1[4];
	q[4] = b2[4]*b1[1]+b2[3]*b1[2]-b2[2]*b1[3]+b2[1]*b1[4];
	
	return;
}




/*
%
%	addEuler121(E1,E2,Q) computes the overall (1-2-1) Euler
%	angle vector corresponding to two successive
%	(1-2-1) rotations E1 and E2.
%
*/
void addEuler121(double *e1, double *e2, double *q) 
{
	double cp1,cp2,sp1,sp2,dum,cp3;
	
	cp1 = cos(e1[2]);
	cp2 = cos(e2[2]);
	sp1 = sin(e1[2]);
	sp2 = sin(e2[2]);
	dum = e1[3]+e2[1];

	q[2] = acos(cp1*cp2-sp1*sp2*cos(dum));
	cp3 = cos(q[2]);
	q[1] = Picheck(e1[1] + atan2(sp1*sp2*sin(dum),cp2-cp3*cp1));
	q[3] = Picheck(e2[3] + atan2(sp1*sp2*sin(dum),cp1-cp3*cp2)); 
	
	return;
}

/*
%
%	addEuler131(E1,E2,Q) computes the overall (1-3-1) Euler
%	angle vector corresponding to two successive
%	(1-3-1) rotations E1 and E2.
%
*/
void addEuler131(double *e1, double *e2, double *q) 
{
	double cp1,cp2,sp1,sp2,dum,cp3;
	
	cp1 = cos(e1[2]);
	cp2 = cos(e2[2]);
	sp1 = sin(e1[2]);
	sp2 = sin(e2[2]);
	dum = e1[3]+e2[1];

	q[2] = acos(cp1*cp2-sp1*sp2*cos(dum));
	cp3 = cos(q[2]);
	q[1] = Picheck(e1[1] + atan2(sp1*sp2*sin(dum),cp2-cp3*cp1));
	q[3] = Picheck(e2[3] + atan2(sp1*sp2*sin(dum),cp1-cp3*cp2)); 
	
	return;
}


/*
%
%	addEuler123(E1,E2,Q) computes the overall (1-2-3) Euler
%	angle vector corresponding to two successive
%	(1-2-3) rotations E1 and E2.
%
*/
void addEuler123(double *e1, double *e2, double *q) 
{
	double C1[4][4],C2[4][4],C[4][4];
	
	Euler1232C(e1,C1);
	Euler1232C(e2,C2);
	MdotM(C2,C1,C);
	C2Euler123(C,q);
	
	return;
}


/*
%
%   addEuler132(E1,E2,Q) computes the overall (1-3-2) Euler
%	angle vector corresponding to two successive
%	(1-3-2) rotations E1 and E2.
%
*/
void addEuler132(double *e1, double *e2, double *q) 
{
	double C1[4][4],C2[4][4],C[4][4];
	
	Euler1322C(e1,C1);
	Euler1322C(e2,C2);
	MdotM(C2,C1,C);
	C2Euler132(C,q);
	
	return;
}






/*
%
%	addEuler212(E1,E2,Q) computes the overall (2-1-2) Euler
%	angle vector corresponding to two successive
%	(2-1-2) rotations E1 and E2.
%
*/
void addEuler212(double *e1, double *e2, double *q) 
{
	double cp1,cp2,sp1,sp2,dum,cp3;
	
	cp1 = cos(e1[2]);
	cp2 = cos(e2[2]);
	sp1 = sin(e1[2]);
	sp2 = sin(e2[2]);
	dum = e1[3]+e2[1];

	q[2] = acos(cp1*cp2-sp1*sp2*cos(dum));
	cp3 = cos(q[2]);
	q[1] = Picheck(e1[1] + atan2(sp1*sp2*sin(dum),cp2-cp3*cp1));
	q[3] = Picheck(e2[3] + atan2(sp1*sp2*sin(dum),cp1-cp3*cp2)); 
	
	return;
}





/*
%
%	addEuler213(E1,E2,Q) computes the overall (2-1-3) Euler
%	angle vector corresponding to two successive
%	(2-1-3) rotations E1 and E2.
%
*/
void addEuler213(double *e1, double *e2, double *q) 
{
	double C1[4][4],C2[4][4],C[4][4];
	
	Euler2132C(e1,C1);
	Euler2132C(e2,C2);
	MdotM(C2,C1,C);
	C2Euler213(C,q);
	
	return;
}




/*
%
%	addEuler231(E1,E2,Q) computes the overall (2-3-1) Euler
%	angle vector corresponding to two successive
%	(2-3-1) rotations E1 and E2.
%
*/
void addEuler231(double *e1, double *e2, double *q) 
{
	double C1[4][4],C2[4][4],C[4][4];
	
	Euler2312C(e1,C1);
	Euler2312C(e2,C2);
	MdotM(C2,C1,C);
	C2Euler231(C,q);
	
	return;
}







/*
%
%	addEuler232(E1,E2,Q) computes the overall (2-3-2) Euler
%	angle vector corresponding to two successive
%	(2-3-2) rotations E1 and E2.
%
*/
void addEuler232(double *e1, double *e2, double *q) 
{
	double cp1,cp2,sp1,sp2,dum,cp3;
	
	cp1 = cos(e1[2]);
	cp2 = cos(e2[2]);
	sp1 = sin(e1[2]);
	sp2 = sin(e2[2]);
	dum = e1[3]+e2[1];

	q[2] = acos(cp1*cp2-sp1*sp2*cos(dum));
	cp3 = cos(q[2]);
	q[1] = Picheck(e1[1] + atan2(sp1*sp2*sin(dum),cp2-cp3*cp1));
	q[3] = Picheck(e2[3] + atan2(sp1*sp2*sin(dum),cp1-cp3*cp2)); 
	
	return;
}





/*
%
%	addEuler312(E1,E2,Q) computes the overall (3-1-2) Euler
%	angle vector corresponding to two successive
%	(3-1-2) rotations E1 and E2.
%
*/
void addEuler312(double *e1, double *e2, double *q) 
{
	double C1[4][4],C2[4][4],C[4][4];
	
	Euler3122C(e1,C1);
	Euler3122C(e2,C2);
	MdotM(C2,C1,C);
	C2Euler312(C,q);
	
	return;
}





/*
%
%	addEuler313(E1,E2,Q) computes the overall (3-1-3) Euler
%	angle vector corresponding to two successive
%	(3-1-3) rotations E1 and E2.
%
*/
void addEuler313(double *e1, double *e2, double *q) 
{
	double cp1,cp2,sp1,sp2,dum,cp3;
	
	cp1 = cos(e1[2]);
	cp2 = cos(e2[2]);
	sp1 = sin(e1[2]);
	sp2 = sin(e2[2]);
	dum = e1[3]+e2[1];

	q[2] = acos(cp1*cp2-sp1*sp2*cos(dum));
	cp3 = cos(q[2]);
	q[1] = Picheck(e1[1] + atan2(sp1*sp2*sin(dum),cp2-cp3*cp1));
	q[3] = Picheck(e2[3] + atan2(sp1*sp2*sin(dum),cp1-cp3*cp2)); 
	
	return;
}






/*
%
%	addEuler321(E1,E2,Q) computes the overall (3-2-1) Euler
%	angle vector corresponding to two successive
%	(3-2-1) rotations E1 and E2.
%
*/
void addEuler321(double *e1, double *e2, double *q) 
{
	double C1[4][4],C2[4][4],C[4][4];
	
	Euler3212C(e1,C1);
	Euler3212C(e2,C2);
	MdotM(C2,C1,C);
	C2Euler321(C,q);
	
	return;
}





/*
%
%	addEuler323(E1,E2,Q) computes the overall (3-2-3) Euler
%	angle vector corresponding to two successive
%	(3-2-3) rotations E1 and E2.
%
*/
void addEuler323(double *e1, double *e2, double *q) 
{
	double cp1,cp2,sp1,sp2,dum,cp3;
	
	cp1 = cos(e1[2]);
	cp2 = cos(e2[2]);
	sp1 = sin(e1[2]);
	sp2 = sin(e2[2]);
	dum = e1[3]+e2[1];

	q[2] = acos(cp1*cp2-sp1*sp2*cos(dum));
	cp3 = cos(q[2]);
	q[1] = Picheck(e1[1] + atan2(sp1*sp2*sin(dum),cp2-cp3*cp1));
	q[3] = Picheck(e2[3] + atan2(sp1*sp2*sin(dum),cp1-cp3*cp2)); 
	
	return;
}




/*
%
%	Q = addGibbs(Q1,Q2) provides the Gibbs vector
%	which corresponds to performing to successive
%	rotations Q1 and Q2.
%
*/
void addGibbs(double *q1, double *q2, double *q) 
{
	double v1[4],v2[4];
	
	cross(q1,q2,v1);
	add(q2,v1,v2);
	add(q1,v2,v1);
	mult(1./(1.-dot(q1,q2)),v1,q);
	
	return;
}




/*
%
%	addMRP(Q1,Q2,Q) provides the MRP vector
%	which corresponds to performing to successive
%	rotations Q1 and Q2.
%
*/
void addMRP(double *q1, double *q2, double *q)
{
	double v1[4], v2[4];
	
	cross(q1,q2,v1);
	mult(2.,v1,v2);
	mult(1-dot(q2,q2),q1,q);
	add(q,v2,q);
	mult(1-dot(q1,q1),q2,v1);
	add(q,v1,q);
	mult(1/(1+dot(q1,q1)*dot(q2,q2)-2*dot(q1,q2)), q, q);
		
	return;
}





/*
%
%	addPRV(Q1,Q2,Q) provides the principal rotation vector
%	which corresponds to performing to successive
%	prinicipal rotations Q1 and Q2.
%
*/
void addPRV(double *qq1, double *qq2, double *q)
{
	double cp1,cp2,sp1,sp2,p,sp, e1[4], e2[4], q1[5], q2[5];
	
	PRV2elem(qq1,q1);
	PRV2elem(qq2,q2);
	cp1 = cos(q1[1]/2.);
	cp2 = cos(q2[1]/2.); 
	sp1 = sin(q1[1]/2.);
	sp2 = sin(q2[1]/2.);
	set(q1[2],q1[3],q1[4],e1);
	set(q2[2],q2[3],q2[4],e2);	
	
	p = 2*acos(cp1*cp2-sp1*sp2*dot(e1,e2));
	sp = sin(p/2.);
	mult(cp1*sp2,e2,q1);
	mult(cp2*sp1,e1,q2);
	add(q1,q2,q);
	cross(e1,e2,q1);
	mult(sp1*sp2,q1,q2);
	add(q,q2,q);
	mult(p/sp,q,q);
	
	return;
}





/*
%
%	BinvEP(Q,B) returns the 3x4 matrix which relates 
%	the derivative of Euler parameter vector Q to the 
%	body angular velocity vector w.  
%	
%		w = 2 [B(Q)]^(-1) dQ/dt
%	
*/
void BinvEP(double *q, double B[4][5]) 
{
  B[1][1] = - q[2];
  B[1][2] = q[1];
  B[1][3] = q[4];
  B[1][4] = -q[3];
  B[2][1] = -q[3];
  B[2][2] = -q[4];
  B[2][3] = q[1];
  B[2][4] = q[2];
  B[3][1] = -q[4];
  B[3][2] = q[3];
  B[3][3] = -q[2];
  B[3][4] = q[1];
	
  return;
}




/*
%
%	BinvEuler121(Q,B) returns the 3x3 matrix which relates 
%	the derivative of the (1-2-1) Euler angle vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	
*/
void BinvEuler121(double *q, double B[4][4]) 
{
  double s2, c2, s3, c3;
  
  s2 = sin(q[2]);
  c2 = cos(q[2]);
  s3 = sin(q[3]);
  c3 = cos(q[3]);

  B[1][1] = c2;
  B[1][2] = 0;
  B[1][3] = 1;
  B[2][1] = s2*s3;
  B[2][2] = c3;
  B[2][3] = 0;
  B[3][1] = s2*c3;
  B[3][2] = -s3;
  B[3][3] = 0;
  
  return;
}





/*
%
%	BinvEuler123(Q,B) returns the 3x3 matrix which relates 
%	the derivative of the (1-2-3) Euler angle vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	
*/
void BinvEuler123(double *q, double B[4][4]) 
{
  double s2, c2, s3, c3;

  s2 = sin(q[2]);
  c2 = cos(q[2]);
  s3 = sin(q[3]);
  c3 = cos(q[3]);

  B[1][1] = c2*c3;
  B[1][2] = s3;
  B[1][3] = 0;
  B[2][1] = -c2*s3;
  B[2][2] = c3;
  B[2][3] = 0;
  B[3][1] = s2;
  B[3][2] = 0;
  B[3][3] = 1;
  
  return;
}





/*
%	BinvEuler131(Q,B) returns the 3x3 matrix which relates 
%	the derivative of the (1-3-1) Euler angle vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	
*/
void BinvEuler131(double *q, double B[4][4]) 
{
  double s2, c2, s3, c3;

  s2 = sin(q[2]);
  c2 = cos(q[2]);
  s3 = sin(q[3]);
  c3 = cos(q[3]);

  B[1][1] = c2;
  B[1][2] = 0;
  B[1][3] = 1;
  B[2][1] = -s2*c3;
  B[2][2] = s3;
  B[2][3] = 0;
  B[3][1] = s2*s3;
  B[3][2] = c3;
  B[3][3] = 0;
  
  return;
}






/*
%
%	BinvEuler132(Q,B) returns the 3x3 matrix which relates 
%	the derivative of the (1-3-2) Euler angle vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	
*/
void BinvEuler132(double *q, double B[4][4]) 
{
  double s2, c2, s3, c3;

  s2 = sin(q[2]);
  c2 = cos(q[2]);
  s3 = sin(q[3]);
  c3 = cos(q[3]);

  B[1][1] = c2*c3;
  B[1][2] = -s3;
  B[1][3] = 0;
  B[2][1] = -s2;
  B[2][2] = 0;
  B[2][3] = 1;
  B[3][1] = c2*s3;
  B[3][2] = c3;
  B[3][3] = 0;
  
  return;
}






/*
%
%	BinvEuler212(Q,B) returns the 3x3 matrix which relates 
%	the derivative of the (2-1-2) Euler angle vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	
*/
void BinvEuler212(double *q, double B[4][4]) 
{
  double s2, c2, s3, c3;

  s2 = sin(q[2]);
  c2 = cos(q[2]);
  s3 = sin(q[3]);
  c3 = cos(q[3]);

  B[1][1] = s2*s3;
  B[1][2] = c3;
  B[1][3] = 0;
  B[2][1] = c2;
  B[2][2] = 0;
  B[2][3] = 1;
  B[3][1] = -s2*c3;
  B[3][2] = s3;
  B[3][3] = 0;
  
  return;
}






/*
%
%	BinvEuler213(Q,B) returns the 3x3 matrix which relates 
%	the derivative of the (2-1-3) Euler angle vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	
*/
void BinvEuler213(double *q, double B[4][4]) 
{
  double s2, c2, s3, c3;

  s2 = sin(q[2]);
  c2 = cos(q[2]);
  s3 = sin(q[3]);
  c3 = cos(q[3]);

  B[1][1] = c2*s3;
  B[1][2] = c3;
  B[1][3] = 0;
  B[2][1] = c2*c3;
  B[2][2] = -s3;
  B[2][3] = 0;
  B[3][1] = -s2;
  B[3][2] = 0;
  B[3][3] = 1;
  
  return;
}






/*
%
%	BinvEuler231(Q,B) returns the 3x3 matrix which relates 
%	the derivative of the (2-3-1) Euler angle vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	
*/
void BinvEuler231(double *q, double B[4][4]) 
{
  double s2, c2, s3, c3;

  s2 = sin(q[2]);
  c2 = cos(q[2]);
  s3 = sin(q[3]);
  c3 = cos(q[3]);

  B[1][1] = s2;
  B[1][2] = 0;
  B[1][3] = 1;
  B[2][1] = c2*c3;
  B[2][2] = s3;
  B[2][3] = 0;
  B[3][1] = -c2*s3;
  B[3][2] = c3;
  B[3][3] = 0;
  
  return;
}







/*
%
%	BinvEuler232(Q,B) returns the 3x3 matrix which relates 
%	the derivative of the (2-3-2) Euler angle vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	
*/
void BinvEuler232(double *q, double B[4][4]) 
{
  double s2, c2, s3, c3;

  s2 = sin(q[2]);
  c2 = cos(q[2]);
  s3 = sin(q[3]);
  c3 = cos(q[3]);

  B[1][1] = s2*c3;
  B[1][2] = -s3;
  B[1][3] = 0;
  B[2][1] = c2;
  B[2][2] = 0;
  B[2][3] = 1;
  B[3][1] = s2*s3;
  B[3][2] = c3;
  B[3][3] = 0;
  
  return;
}








/*
%
%	BinvEuler323(Q,B) returns the 3x3 matrix which relates 
%	the derivative of the (3-2-3) Euler angle vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	
*/
void BinvEuler323(double *q, double B[4][4]) 
{
  double s2, c2, s3, c3;

  s2 = sin(q[2]);
  c2 = cos(q[2]);
  s3 = sin(q[3]);
  c3 = cos(q[3]);

  B[1][1] = -s2*c3;
  B[1][2] = s3;
  B[1][3] = 0;
  B[2][1] = s2*s3;
  B[2][2] = c3;
  B[2][3] = 0;
  B[3][1] = c2;
  B[3][2] = 0;
  B[3][3] = 1;
  
  return;
}








/*
%
%	BinvEuler313(Q,B) returns the 3x3 matrix which relates 
%	the derivative of the (3-1-3) Euler angle vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	
*/
void BinvEuler313(double *q, double B[4][4]) 
{
  double s2, c2, s3, c3;

  s2 = sin(q[2]);
  c2 = cos(q[2]);
  s3 = sin(q[3]);
  c3 = cos(q[3]);

  B[1][1] = s2*s3;
  B[1][2] = c3;
  B[1][3] = 0;
  B[2][1] = s2*c3;
  B[2][2] = -s3;
  B[2][3] = 0;
  B[3][1] = c2;
  B[3][2] = 0;
  B[3][3] = 1;
  
  return;
}








/*
%
%	BinvEuler321(Q,B) returns the 3x3 matrix which relates 
%	the derivative of the (3-2-1) Euler angle vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	
*/
void BinvEuler321(double *q, double B[4][4]) 
{
  double s2, c2, s3, c3;

  s2 = sin(q[2]);
  c2 = cos(q[2]);
  s3 = sin(q[3]);
  c3 = cos(q[3]);

  B[1][1] = -s2;
  B[1][2] = 0;
  B[1][3] = 1;
  B[2][1] = c2*s3;
  B[2][2] = c3;
  B[2][3] = 0;
  B[3][1] = c2*c3;
  B[3][2] = -s3;
  B[3][3] = 0;
  
  return;
}






/*
%
%	BinvEuler312(Q) returns the 3x3 matrix which relates 
%	the derivative of the (3-2-3) Euler angle vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	
*/
void BinvEuler312(double *q, double B[4][4]) 
{
  double s2, c2, s3, c3;

  s2 = sin(q[2]);
  c2 = cos(q[2]);
  s3 = sin(q[3]);
  c3 = cos(q[3]);

  B[1][1] = -c2*s3;
  B[1][2] = c3;
  B[1][3] = 0;
  B[2][1] = s2;
  B[2][2] = 0;
  B[2][3] = 1;
  B[3][1] = c2*c3;
  B[3][2] = s3;
  B[3][3] = 0;
  
  return;
}





/*
%
%	BinvGibbs(Q,B) returns the 3x3 matrix which relates 
%	the derivative of Gibbs vector Q to the 
%	body angular velocity vector w.  
%	
%		w = 2 [B(Q)]^(-1) dQ/dt
%
*/
void BinvGibbs(double *q, double B[4][4])	
{
  B[1][1] = 1;
  B[1][2] = q[3];
  B[1][3] = -q[2];
  B[2][1] = -q[3];
  B[2][2] = 1;
  B[2][3] = q[1];
  B[3][1] = q[2];
  B[3][2] = -q[1];
  B[3][3] = 1;
  Mmult(1./(1+dot(q,q)),B,B);
  
  return;
}





/*
%
%	BinvMRP(Q,B) returns the 3x3 matrix which relates 
%	the derivative of MRP vector Q to the 
%	body angular velocity vector w.  
%	
%		w = 4 [B(Q)]^(-1) dQ/dt
%	
*/
void BinvMRP(double *q, double B[4][4])
{
  double s2;
  
  s2 = dot(q,q);
  B[1][1] = 1-s2+2*q[1]*q[1];
  B[1][2] = 2*(q[1]*q[2]+q[3]);
  B[1][3] = 2*(q[1]*q[3]-q[2]);
  B[2][1] = 2*(q[2]*q[1]-q[3]);
  B[2][2] = 1-s2+2*q[2]*q[2];
  B[2][3] = 2*(q[2]*q[3]+q[1]);
  B[3][1] = 2*(q[3]*q[1]+q[2]);
  B[3][2] = 2*(q[3]*q[2]-q[1]);
  B[3][3] = 1-s2+2*q[3]*q[3];
  Mmult(1./(1+s2)/(1+s2),B,B);

  return;
}





/*
%
%	BinvPRV(Q,B) returns the 3x3 matrix which relates 
%	the derivative of principal rotation vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	
*/
void BinvPRV(double *q, double B[4][4])
{
  double p, c1, c2;
  
  p = sqrt(dot(q,q));
  c1 = (1-cos(p))/p/p;
  c2 = (p-sin(p))/p/p/p;

  B[1][1] = 1-c2*(q[2]*q[2]+q[3]*q[3]);
  B[1][2] = c1*q[3] + c2*q[1]*q[2];
  B[1][3] = -c1*q[2] + c2*q[1]*q[3];
  B[2][1] = -c1*q[3] + c2*q[1]*q[2];
  B[2][2] = 1 - c2*(q[1]*q[1]+ q[3]*q[3]);
  B[2][3] = c1*q[1] + c2*q[2]*q[3];
  B[3][1] = c1*q[2] + c2*q[3]*q[1];
  B[3][2] = -c1*q[1] + c2*q[3]*q[2];
  B[3][3] = 1 - c2*(q[1]*q[1]+q[2]*q[2]);
  
  return;
}





/*
%
%	BmatEP(Q,B) returns the 4x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	Euler parameter vector Q.  
%	
%		dQ/dt = 1/2 [B(Q)] w
%	
*/
void BmatEP(double *q, double B[5][4])
{
  
  B[1][1] = -q[2];
  B[1][2] = -q[3];
  B[1][3] = -q[4];
  B[2][1] = q[1];
  B[2][2] = -q[4];
  B[2][3] = q[3];
  B[3][1] = q[4];
  B[3][2] = q[1];
  B[3][3] = -q[2];
  B[4][1] = -q[3];
  B[4][2] = q[2];
  B[4][3] = q[1];
  
  return;
}



/*
%
%	BmatEuler121(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	(1-2-1) Euler angle vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	
*/
void BmatEuler121(double *q, double B[4][4])
{
    double s2, c2, s3, c3;
    
    s2 = sin(q[2]);
    c2 = cos(q[2]);
    s3 = sin(q[3]);
    c3 = cos(q[3]);

    B[1][1] = 0;
    B[1][2] = s3;
    B[1][3] = c3;
    B[2][1] = 0;
    B[2][2] = s2*c3;
    B[2][3] = -s2*s3;
    B[3][1] = s2;
    B[3][2] = -c2*s3;
    B[3][3] = -c2*c3;
    Mmult(1./s2,B,B);
    
    return;
}





/*
%
%	BmatEuler131(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	(1-3-1) Euler angle vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	
*/
void BmatEuler131(double *q, double B[4][4])
{
    double s2, c2, s3, c3;
    
    s2 = sin(q[2]);
    c2 = cos(q[2]);
    s3 = sin(q[3]);
    c3 = cos(q[3]);

    B[1][1] = 0;
    B[1][2] = -c3;
    B[1][3] = s3;
    B[2][1] = 0;
    B[2][2] = s2*s3;
    B[2][3] = s2*c3;
    B[3][1] = s2;
    B[3][2] = c2*c3;
    B[3][3] = -c2*s3;
    Mmult(1./s2,B,B);
    
    return;
}




/*
%
%	BmatEuler123(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	(1-2-3) Euler angle vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	
*/
void BmatEuler123(double *q, double B[4][4])
{
    double s2, c2, s3, c3;
    
    s2 = sin(q[2]);
    c2 = cos(q[2]);
    s3 = sin(q[3]);
    c3 = cos(q[3]);

    B[1][1] = c3;
    B[1][2] = -s3;
    B[1][3] = 0;
    B[2][1] = c2*s3;
    B[2][2] = c2*c3;
    B[2][3] = 0;
    B[3][1] = -s2*c3;
    B[3][2] = s2*s3;
    B[3][3] = c2;
    Mmult(1./c2,B,B);
    
    return;
}







/*
%
%	BmatEuler132(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	(1-3-2) Euler angle vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	
*/
void BmatEuler132(double *q, double B[4][4])
{
    double s2, c2, s3, c3;
    
    s2 = sin(q[2]);
    c2 = cos(q[2]);
    s3 = sin(q[3]);
    c3 = cos(q[3]);

    B[1][1] = c3;
    B[1][2] = 0;
    B[1][3] = s3;
    B[2][1] = -c2*s3;
    B[2][2] = 0;
    B[2][3] = c2*c3;
    B[3][1] = s2*c3;
    B[3][2] = c2;
    B[3][3] = s2*s3;
    Mmult(1./c2,B,B);
    
    return;
}






/*
%
%	BmatEuler212(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	(2-1-2) Euler angle vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	
*/
void BmatEuler212(double *q, double B[4][4])
{
    double s2, c2, s3, c3;
    
    s2 = sin(q[2]);
    c2 = cos(q[2]);
    s3 = sin(q[3]);
    c3 = cos(q[3]);

    B[1][1] = s3;
    B[1][2] = 0;
    B[1][3] = -c3;
    B[2][1] = s2*c3;
    B[2][2] = 0;
    B[2][3] = s2*s3;
    B[3][1] = -c2*s3;
    B[3][2] = s2;
    B[3][3] = c2*c3;
    Mmult(1./s2,B,B);
    
    return;
}






/*
%
%	BmatEuler213(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	(2-1-3) Euler angle vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	
*/
void BmatEuler213(double *q, double B[4][4])
{
    double s2, c2, s3, c3;
    
    s2 = sin(q[2]);
    c2 = cos(q[2]);
    s3 = sin(q[3]);
    c3 = cos(q[3]);

    B[1][1] = s3;
    B[1][2] = c3;
    B[1][3] = 0;
    B[2][1] = c2*c3;
    B[2][2] = -c2*s3;
    B[2][3] = 0;
    B[3][1] = s2*s3;
    B[3][2] = s2*c3;
    B[3][3] = c2;
    Mmult(1./c2,B,B);
    
    return;
}






/*
%
%	BmatEuler231(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	(2-3-1) Euler angle vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	
*/
void BmatEuler231(double *q, double B[4][4])
{
    double s2, c2, s3, c3;
    
    s2 = sin(q[2]);
    c2 = cos(q[2]);
    s3 = sin(q[3]);
    c3 = cos(q[3]);

    B[1][1] = 0;
    B[1][2] = c3;
    B[1][3] = -s3;
    B[2][1] = 0;
    B[2][2] = c2*s3;
    B[2][3] = c2*c3;
    B[3][1] = c2;
    B[3][2] = -s2*c3;
    B[3][3] = s2*s3;
    Mmult(1./c2,B,B);
    
    return;
}







/*
%
%	B = BmatEuler232(Q) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	(2-3-2) Euler angle vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	
*/
void BmatEuler232(double *q, double B[4][4])
{
    double s2, c2, s3, c3;
    
    s2 = sin(q[2]);
    c2 = cos(q[2]);
    s3 = sin(q[3]);
    c3 = cos(q[3]);

    B[1][1] = c3;
    B[1][2] = 0;
    B[1][3] = s3;
    B[2][1] = -s2*s3;
    B[2][2] = 0;
    B[2][3] = s2*c3;
    B[3][1] = -c2*c3;
    B[3][2] = s2;
    B[3][3] = -c2*s3;
    Mmult(1./s2,B,B);
    
    return;
}







/*
%
%	BmatEuler312(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	(3-1-2) Euler angle vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	
*/
void BmatEuler312(double *q, double B[4][4])
{
    double s2, c2, s3, c3;
    
    s2 = sin(q[2]);
    c2 = cos(q[2]);
    s3 = sin(q[3]);
    c3 = cos(q[3]);

    B[1][1] = -s3;
    B[1][2] = 0;
    B[1][3] = c3;
    B[2][1] = c2*c3;
    B[2][2] = 0;
    B[2][3] = c2*s3;
    B[3][1] = s2*s3;
    B[3][2] = c2;
    B[3][3] = -s2*c3;
    Mmult(1./c2,B,B);
    
    return;
}







/*
%
%	BmatEuler313(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	(3-1-3) Euler angle vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	
*/
void BmatEuler313(double *q, double B[4][4])
{
    double s2, c2, s3, c3;
    
    s2 = sin(q[2]);
    c2 = cos(q[2]);
    s3 = sin(q[3]);
    c3 = cos(q[3]);

    B[1][1] = s3;
    B[1][2] = c3;
    B[1][3] = 0;
    B[2][1] = c3*s2;
    B[2][2] = -s3*s2;
    B[2][3] = 0;
    B[3][1] = -s3*c2;
    B[3][2] = -c3*c2;
    B[3][3] = s2;
    Mmult(1./s2,B,B);
    
    return;
}






/*
%
%	BmatEuler321(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	(3-2-1) Euler angle vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	
*/
void BmatEuler321(double *q, double B[4][4])
{
    double s2, c2, s3, c3;
    
    s2 = sin(q[2]);
    c2 = cos(q[2]);
    s3 = sin(q[3]);
    c3 = cos(q[3]);

    B[1][1] = 0;
    B[1][2] = s3;
    B[1][3] = c3;
    B[2][1] = 0;
    B[2][2] = c2*c3;
    B[2][3] = -c2*s3;
    B[3][1] = c2;
    B[3][2] = s2*s3;
    B[3][3] = s2*c3;
    Mmult(1./c2,B,B);
    
    return;
}







/*
%
%	BmatEuler323(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	(3-2-3) Euler angle vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	
*/
void BmatEuler323(double *q, double B[4][4])
{
    double s2, c2, s3, c3;
    
    s2 = sin(q[2]);
    c2 = cos(q[2]);
    s3 = sin(q[3]);
    c3 = cos(q[3]);

    B[1][1] = -c3;
    B[1][2] = s3;
    B[1][3] = 0;
    B[2][1] = s2*s3;
    B[2][2] = s2*c3;
    B[2][3] = 0;
    B[3][1] = c2*c3;
    B[3][2] = -c2*s3;
    B[3][3] = s2;
    Mmult(1./s2,B,B);
    
    return;
}







/*
%
%	BmatGibbs(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	Gibbs vector Q.  
%	
%		dQ/dt = 1/2 [B(Q)] w
%	
*/
void BmatGibbs(double *q, double B[4][4])
{

    B[1][1] = 1+q[1]*q[1];
    B[1][2] = q[1]*q[2]-q[3];
    B[1][3] = q[1]*q[3]+q[2];
    B[2][1] = q[2]*q[1]+q[3];
    B[2][2] = 1+q[2]*q[2];
    B[2][3] = q[2]*q[3]-q[1];
    B[3][1] = q[3]*q[1]-q[2];
    B[3][2] = q[3]*q[2]+q[1];
    B[3][3] = 1+q[3]*q[3];
    
    return;
}





/*
%
%	BmatMRP(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	MRP vector Q.  
%	
%		dQ/dt = 1/4 [B(Q)] w
%	
*/
void BmatMRP(double *q, double B[4][4])
{
    double s2;
    
    s2 = dot(q,q);
    B[1][1] = 1-s2+2*q[1]*q[1];
    B[1][2] = 2*(q[1]*q[2]-q[3]);
    B[1][3] = 2*(q[1]*q[3]+q[2]);
    B[2][1] = 2*(q[2]*q[1]+q[3]);
    B[2][2] = 1-s2+2*q[2]*q[2];
    B[2][3] = 2*(q[2]*q[3]-q[1]);
    B[3][1] = 2*(q[3]*q[1]-q[2]);
    B[3][2] = 2*(q[3]*q[2]+q[1]);
    B[3][3] = 1-s2+2*q[3]*q[3];
    
    return;
}





/*
%
%	BmatPRV(Q,B) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	principal rotation vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	
*/
void BmatPRV(double *q, double B[4][4])
{
    double p,c;
    p = norm(q);
    c = 1./p/p*(1.-p/2./tan(p/2.));
    B[1][1] = 1- c*(q[2]*q[2]+q[3]*q[3]);
    B[1][2] = -q[3]/2 + c*(q[1]*q[2]);
    B[1][3] = q[2]/2 + c*(q[1]*q[3]);
    B[2][1] = q[3]/2 + c*(q[1]*q[2]);
    B[2][2] = 1 - c*(q[1]*q[1]+q[3]*q[3]);
    B[2][3] = -q[1]/2 + c*(q[2]*q[3]);
    B[3][1] = -q[2]/2 + c*(q[1]*q[3]);
    B[3][2] = q[1]/2 + c*(q[2]*q[3]);
    B[3][3] = 1-c*(q[1]*q[1]+q[2]*q[2]);
    
    return;
}





/*
%
%	C2EP(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding 4x1 Euler parameter vector Q,
%	where the first component of Q is the non-dimensional
%	Euler parameter Beta_0 >= 0. Transformation is done
%	using the Stanley method.i
%	
*/
void    C2EP(double C[4][4], double b[5])
{
    double tr, b2[5],max;
    int i,j;
    
    tr = C[1][1]+C[2][2]+C[3][3];
    b2[1] = (1+tr)/4.;
    b2[2] = (1+2*C[1][1]-tr)/4.;
    b2[3] = (1+2*C[2][2]-tr)/4.;
    b2[4] = (1+2*C[3][3]-tr)/4.;

    i = 1;
    max = b2[1];
    for (j=2;j<=4;j++) {
        if (b2[j] > max) {
            i = j;
            max = b2[j];
        }
    }
    
    switch (i) {
        case 1:
            b[1] = sqrt(b2[1]);
            b[2] = (C[2][3]-C[3][2])/4/b[1];
            b[3] = (C[3][1]-C[1][3])/4/b[1];
            b[4] = (C[1][2]-C[2][1])/4/b[1];
            break;
        case 2:
            b[2] = sqrt(b2[2]);
            b[1] = (C[2][3]-C[3][2])/4/b[2];
            if (b[1]<0) {
                b[2] = -b[2];
                b[1] = -b[1];
            }
            b[3] = (C[1][2]+C[2][1])/4/b[2];
            b[4] = (C[3][1]+C[1][3])/4/b[2];
            break;
        case 3:
            b[3] = sqrt(b2[3]);
            b[1] = (C[3][1]-C[1][3])/4/b[3];
            if (b[1]<0) {
                b[3] = -b[3];
                b[1] = -b[1];
            }
            b[2] = (C[1][2]+C[2][1])/4/b[3];
            b[4] = (C[2][3]+C[3][2])/4/b[3];
            break;
        case 4:
            b[4] = sqrt(b2[4]);
            b[1] = (C[1][2]-C[2][1])/4/b[4];
            if (b[1]<0) {
                b[4] = -b[4];
                b[1] = -b[1];
            }
            b[2] = (C[3][1]+C[1][3])/4/b[4];
            b[3] = (C[2][3]+C[3][2])/4/b[4];
            break;
    }
    
    return;
}



/*
%
%	C2Euler121(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding (1-2-1) Euler angle set.
%	
*/
void    C2Euler121(double C[4][4], double *q)
{
    q[1] = atan2(C[1][2],-C[1][3]);
    q[2] = acos(C[1][1]);
    q[3]= atan2(C[2][1],C[3][1]);
    
    return;
}





/*
%
%	C2Euler123(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding (1-2-3) Euler angle set.
%	
*/
void    C2Euler123(double C[4][4], double *q)
{
    q[1] = atan2(-C[3][2],C[3][3]);
    q[2] = asin(C[3][1]);
    q[3]= atan2(-C[2][1],C[1][1]);
    
    return;
}





/*
%
%	C2Euler131(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding (1-3-1) Euler angle set.
%	
*/
void    C2Euler131(double C[4][4], double *q)
{
    q[1] = atan2(C[1][3],C[1][2]);
    q[2] = acos(C[1][1]);
    q[3]= atan2(C[3][1],-C[2][1]);
    
    return;
}




/*
%
%	C2Euler132(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding (1-3-2) Euler angle set.
%	
*/
void    C2Euler132(double C[4][4], double *q)
{
    q[1] = atan2(C[2][3],C[2][2]);
    q[2] = asin(-C[2][1]);
    q[3] = atan2(C[3][1],C[1][1]);
    
    return;
}





/*
%
%	C2Euler212(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding (2-1-2) Euler angle set.
%	
*/
void    C2Euler212(double C[4][4], double *q)
{
    q[1] = atan2(C[2][1],C[2][3]);
    q[2] = acos(C[2][2]);
    q[3] = atan2(C[1][2],-C[3][2]);
    
    return;
}




/*
%
%	C2Euler213(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding (2-1-3) Euler angle set.
%	
*/
void    C2Euler213(double C[4][4], double *q)
{
    q[1] = atan2(C[3][1],C[3][3]);
    q[2] = asin(-C[3][2]);
    q[3]= atan2(C[1][2],C[2][2]);

    return;
}




/*
%
%	C2Euler231(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding (2-3-1) Euler angle set.
%	
*/
void    C2Euler231(double C[4][4], double *q)
{
    q[1] = atan2(-C[1][3],C[1][1]);
    q[2] = asin(C[1][2]);
    q[3]= atan2(-C[3][2],C[2][2]);
    
    return;
}





/*
%
%	C2Euler232(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding (2-3-2) Euler angle set.
%	
*/
void    C2Euler232(double C[4][4], double *q)
{
    q[1] = atan2(C[2][3],-C[2][1]);
    q[2] = acos(C[2][2]);
    q[3]= atan2(C[3][2],C[1][2]);
    
    return;
}




/*
%
%	C2Euler312(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding (3-1-2) Euler angle set.
%	
*/
void C2Euler312(double C[4][4], double *q)
{
    q[1] = atan2(-C[2][1],C[2][2]);
    q[2] = asin(C[2][3]);
    q[3]= atan2(-C[1][3],C[3][3]);
    
    return;
}





/*
%
%	C2Euler313(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding (3-1-3) Euler angle set.
%	
*/
void    C2Euler313(double C[4][4], double *q)
{
    q[1] = atan2(C[3][1],-C[3][2]);
    q[2] = acos(C[3][3]);
    q[3]= atan2(C[1][3],C[2][3]);

    return;
}




/*
%
%	C2Euler321(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding (3-2-1) Euler angle set.
%	
*/
void    C2Euler321(double C[4][4], double *q)
{
    q[1] = atan2(C[1][2],C[1][1]);
    q[2] = asin(-C[1][3]);
    q[3]= atan2(C[2][3],C[3][3]);
    
    return;
}





/*
%
%	C2Euler323(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding (3-2-3) Euler angle set.
%	
*/
void    C2Euler323(double C[4][4], double *q)
{
    q[1] = atan2(C[3][2],C[3][1]);
    q[2] = acos(C[3][3]);
    q[3]= atan2(C[2][3],-C[1][3]);
    
    return;
}




/*
%
%	C2Gibbs(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding 3x1 Gibbs vector Q.
%	
*/
void    C2Gibbs(double C[4][4], double *q)
{
    double b[5];
    
    C2EP(C,b);

    q[1] = b[2]/b[1];
    q[2] = b[3]/b[1];
    q[3] = b[4]/b[1];
    
    return;
}






/*
%
%	C2MRP(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding 3x1 MRP vector Q where the 
%	MRP vector is chosen such that |Q| <= 1.
%
*/
void    C2MRP(double C[4][4], double *q)
{
    double b[5];
    
    C2EP(C,b);

    q[1] = b[2]/(1+b[1]);
    q[2] = b[3]/(1+b[1]);
    q[3] = b[4]/(1+b[1]);
    
    return;
}




/*
%
%	C2PRV(C,Q) translates the 3x3 direction cosine matrix
%	C into the corresponding 3x1 principal rotation vector Q,
%	where the first component of Q is the principal rotation angle
%	phi (0<= phi <= Pi)
%	
*/
void    C2PRV(double C[4][4], double *q)
{
    double cp, p, sp;

    cp = (C[1][1] + C[2][2] + C[3][3]-1)/2;
    p = acos(cp);
    sp = p/2./sin(p);
    q[1] = (C[2][3]-C[3][2])*sp;
    q[2] = (C[3][1]-C[1][3])*sp;
    q[3] = (C[1][2]-C[2][1])*sp;
    
    return;
}






/*
%
%	dEP(Q,W,dq) returns the Euler parameter derivative
%	for a given Euler parameter vector Q and body
%	angular velocity vector w.
%
%	dQ/dt = 1/2 [B(Q)] w
%
*/
void    dEP(double *q, double *w, double *dq)
{
    double B[5][4];
	int i,j;
    
    BmatEP(q,B);
    Mdot(B,w,dq);
	for (i=1;i<=4;i++) {
		dq[i] = 0.;
		for (j=1;j<=3;j++) {
			dq[i] += B[i][j]*w[j];
		}
	}
    mult(.5,dq,dq);
    dq[4] = 0.5*dq[4];
	
    return;
}




/*
%
%	dEuler121(Q,W,dq) returns the (1-2-1) Euler angle derivative
%	vector for a given (1-2-1) Euler angle vector Q and body
%	angular velocity vector w.
%
%	dQ/dt =  [B(Q)] w
%
*/
void    dEuler121(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatEuler121(q,B);
    Mdot(B,w,dq);
    
    return;
}





/*
%
%	dEuler123(Q,W,dq) returns the (1-2-3) Euler angle derivative
%	vector for a given (1-2-3) Euler angle vector Q and body
%	angular velocity vector w.
%
%	dQ/dt =  [B(Q)] w
%
*/
void    dEuler123(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatEuler123(q,B);
    Mdot(B,w,dq);
    
    return;
}






/*
%
%	dEuler131(Q,W,dq) returns the (1-3-1) Euler angle derivative
%	vector for a given (1-3-1) Euler angle vector Q and body
%	angular velocity vector w.
%
%	dQ/dt =  [B(Q)] w
%
*/
void    dEuler131(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatEuler131(q,B);
    Mdot(B,w,dq);
    
    return;
}






/*
%
%	dEuler132(Q,W,dq) returns the (1-3-2) Euler angle derivative
%	vector for a given (1-3-2) Euler angle vector Q and body
%	angular velocity vector w.
%
%	dQ/dt =  [B(Q)] w
%
*/
void    dEuler132(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatEuler132(q,B);
    Mdot(B,w,dq);
    
    return;
}





/*
%
%	dEuler212(Q,W,dq) returns the (2-1-2) Euler angle derivative
%	vector for a given (2-1-2) Euler angle vector Q and body
%	angular velocity vector w.
%
%	dQ/dt =  [B(Q)] w
%
*/
void    dEuler212(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatEuler212(q,B);
    Mdot(B,w,dq);
    
    return;
}





/*
%
%	dEuler213(Q,W,dq) returns the (2-1-3) Euler angle derivative
%	vector for a given (2-1-3) Euler angle vector Q and body
%	angular velocity vector w.
%
%	dQ/dt =  [B(Q)] w
%
*/
void    dEuler213(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatEuler213(q,B);
    Mdot(B,w,dq);
    
    return;
}





/*
%
%	dEuler231(Q,W,dq) returns the (2-3-1) Euler angle derivative
%	vector for a given (2-3-1) Euler angle vector Q and body
%	angular velocity vector w.
%
%	dQ/dt =  [B(Q)] w
%
*/
void    dEuler231(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatEuler231(q,B);
    Mdot(B,w,dq);
    
    return;
}





/*
%
%	dEuler232(Q,W,dq) returns the (2-3-2) Euler angle derivative
%	vector for a given (2-3-2) Euler angle vector Q and body
%	angular velocity vector w.
%
%	dQ/dt =  [B(Q)] w
%
*/
void    dEuler232(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatEuler232(q,B);
    Mdot(B,w,dq);
    
    return;
}







/*
%
%	dEuler312(Q,W,dq) returns the (3-1-2) Euler angle derivative
%	vector for a given (3-1-2) Euler angle vector Q and body
%	angular velocity vector w.
%
%	dQ/dt =  [B(Q)] w
%
*/
void    dEuler312(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatEuler312(q,B);
    Mdot(B,w,dq);
    
    return;
}






/*
%
%	dEuler313(Q,W,dq) returns the (3-1-3) Euler angle derivative
%	vector for a given (3-1-3) Euler angle vector Q and body
%	angular velocity vector w.
%
%	dQ/dt =  [B(Q)] w
%
*/
void    dEuler313(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatEuler313(q,B);
    Mdot(B,w,dq);
    
    return;
}





/*
%
%	dEuler321(Q,W,dq) returns the (3-2-1) Euler angle derivative
%	vector for a given (3-2-1) Euler angle vector Q and body
%	angular velocity vector w.
%
%	dQ/dt =  [B(Q)] w
%
*/
void    dEuler321(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatEuler321(q,B);
    Mdot(B,w,dq);
    
    return;
}






/*
%
%	dEuler323(Q,W,dq) returns the (3-2-3) Euler angle derivative
%	vector for a given (3-2-3) Euler angle vector Q and body
%	angular velocity vector w.
%
%	dQ/dt =  [B(Q)] w
%
*/
void    dEuler323(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatEuler323(q,B);
    Mdot(B,w,dq);
    
    return;
}




/*
%
%	dGibbs(Q,W,dq) returns the Gibbs derivative
%	for a given Gibbs vector Q and body
%	angular velocity vector w.
%
%	dQ/dt = 1/2 [B(Q)] w
%
*/
void    dGibbs(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatGibbs(q,B);
    Mdot(B,w,dq);
    mult(0.5,dq,dq);
    
    return;
}


/*
%
%	dMRP(Q,W,dq) returns the MRP derivative
%	for a given MRP vector Q and body
%	angular velocity vector w.
%
%	dQ/dt = 1/4 [B(Q)] w
%
*/
void    dMRP(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatMRP(q,B);
    Mdot(B,w,dq);
    mult(0.25,dq,dq);
    
    return;
}




/*
%
%	dPRV(Q,W,dq) returns the PRV derivative
%	for a given PRV vector Q and body
%	angular velocity vector w.
%
%	dQ/dt =  [B(Q)] w
%
*/
void    dPRV(double *q, double *w, double *dq)
{
    double B[4][4];
    
    BmatPRV(q,B);
    Mdot(B,w,dq);
    
    return;
}





/*
%	
%	elem2PRV(R,Q) translates a prinicpal rotation 
%	element set R into the corresponding principal 
%	rotation vector Q.
%
*/
void	elem2PRV(double *r, double *q)
{
  q[1] = r[2]*r[1];
  q[2] = r[3]*r[1];
  q[3] = r[4]*r[1];
  
  return;
}




/*
%
%	EP2C(Q,C) returns the direction cosine 
%	matrix in terms of the 4x1 Euler parameter vector
%	Q.  The first element is the non-dimensional Euler
%	parameter, while the remain three elements form 
%	the Eulerparameter vector.
%
*/
void	EP2C(double *q, double C[4][4])
{
  double q0, q1, q2, q3;
  
  q0 = q[1];
  q1 = q[2];
  q2 = q[3];
  q3 = q[4];

  C[1][1] = q0*q0+q1*q1-q2*q2-q3*q3;
  C[1][2] = 2*(q1*q2+q0*q3);
  C[1][3] = 2*(q1*q3-q0*q2);
  C[2][1] = 2*(q1*q2-q0*q3);
  C[2][2] = q0*q0-q1*q1+q2*q2-q3*q3;
  C[2][3] = 2*(q2*q3+q0*q1);
  C[3][1] = 2*(q1*q3 + q0*q2);
  C[3][2] = 2*(q2*q3-q0*q1);
  C[3][3] = q0*q0-q1*q1-q2*q2+q3*q3;
  
  return;
}





/*
%
%	EP2Euler121(Q,E) translates the Euler parameter
%	vector Q into the corresponding (1-2-1) Euler angle
%	vector E.
%
*/
void	EP2Euler121(double *q, double *e)
{
  double t1, t2;
  
  t1 = atan2(q[4],q[3]);
  t2 = atan2(q[2],q[1]);

  e[1] = t1+t2;
  e[2] = 2*acos(sqrt(q[1]*q[1]+q[2]*q[2]));
  e[3] = t2-t1;
  
  return;
}




/*
%
%	EP2Euler123(Q,E) translates the Euler parameter vector
%	Q into the corresponding (1-2-3) Euler angle set.
%	
*/
void	EP2Euler123(double *q, double *e)
{
  double q0, q1, q2, q3;
  
  q0 = q[1];
  q1 = q[2];
  q2 = q[3];
  q3 = q[4];

  e[1] = atan2(-2*(q2*q3-q0*q1),q0*q0-q1*q1-q2*q2+q3*q3);
  e[2] = asin(2*(q1*q3 + q0*q2));
  e[3] = atan2(-2*(q1*q2-q0*q3),q0*q0+q1*q1-q2*q2-q3*q3);
  
  return;
}





/*
%
%	EP2Euler131(Q,E) translates the Euler parameter
%	vector Q into the corresponding (1-3-1) Euler angle
%	vector E.
%
*/
void	EP2Euler131(double *q, double *e)
{
  double t1, t2;

  t1 = atan2(q[3],q[4]);
  t2 = atan2(q[2],q[1]);

  e[1] = t2-t1;
  e[2] = 2*acos(sqrt(q[1]*q[1]+q[2]*q[2]));
  e[3] = t2+t1;
  
  return;
}






/*
%
%	EP2Euler132(Q,E) translates the Euler parameter vector
%	Q into the corresponding (1-3-2) Euler angle set.
%	
*/
void	EP2Euler132(double *q, double *e)
{
  double q0, q1, q2, q3;
  
  q0 = q[1];
  q1 = q[2];
  q2 = q[3];
  q3 = q[4];

  e[1] = atan2(2*(q2*q3+q0*q1),q0*q0-q1*q1+q2*q2-q3*q3);
  e[2] = asin(-2*(q1*q2-q0*q3));
  e[3]= atan2(2*(q1*q3 + q0*q2),q0*q0+q1*q1-q2*q2-q3*q3);
  
  return;
}




/*
%
%	EP2Euler212(Q,E) translates the Euler parameter
%	vector Q into the corresponding (2-1-2) Euler angle
%	vector E.
%
*/
void	EP2Euler212(double *q, double *e)
{
  double t1,t2;

  t1 = atan2(q[4],q[2]);
  t2 = atan2(q[3],q[1]);

  e[1] = t2-t1;
  e[2] = 2*acos(sqrt(q[1]*q[1]+q[3]*q[3]));
  e[3] = t2+t1;
  
  return;
}





/*
%
%	EP2Euler213(Q,E) translates the Euler parameter vector
%	Q into the corresponding (2-1-3) Euler angle set.
%	
*/
void	EP2Euler213(double *q, double *e)
{
  double q0, q1, q2, q3;
  
  q0 = q[1];
  q1 = q[2];
  q2 = q[3];
  q3 = q[4];

  e[1] = atan2(2*(q1*q3 + q0*q2),q0*q0-q1*q1-q2*q2+q3*q3);
  e[2] = asin(-2*(q2*q3-q0*q1));
  e[3]= atan2(2*(q1*q2+q0*q3),q0*q0-q1*q1+q2*q2-q3*q3);
  
  return;
}





/*
%
%	EP2Euler231(Q,E) translates the Euler parameter vector
%	Q into the corresponding (2-3-1) Euler angle set.
%	
*/
void	EP2Euler231(double *q, double *e)
{
  double q0, q1, q2, q3;
  
  q0 = q[1];
  q1 = q[2];
  q2 = q[3];
  q3 = q[4];

  e[1] = atan2(-2*(q1*q3-q0*q2), q0*q0+q1*q1-q2*q2-q3*q3);
  e[2] = asin(2*(q1*q2+q0*q3));
  e[3]= atan2(-2*(q2*q3-q0*q1),q0*q0-q1*q1+q2*q2-q3*q3);
  
  return;
}






/*
%
%	EP2Euler232(Q,E) translates the Euler parameter
%	vector Q into the corresponding (2-3-2) Euler angle
%	vector E.
%
*/
void	EP2Euler232(double *q, double *e)
{
  double t1, t2;
  
  t1 = atan2(q[2],q[4]);
  t2 = atan2(q[3],q[1]);

  e[1] = t1+t2;
  e[2] = 2*acos(sqrt(q[1]*q[1]+q[3]*q[3]));
  e[3] = t2-t1;
  
  return;
}






/*
%
%	EP2Euler312(Q,E) translates the Euler parameter vector
%	Q into the corresponding (3-1-2) Euler angle set.
%	
*/
void	EP2Euler312(double *q, double *e)
{
  double q0, q1, q2, q3;
  
  q0 = q[1];
  q1 = q[2];
  q2 = q[3];
  q3 = q[4];

  e[1] = atan2(-2*(q1*q2-q0*q3),q0*q0-q1*q1+q2*q2-q3*q3);
  e[2] = asin(2*(q2*q3+q0*q1));
  e[3]= atan2(-2*(q1*q3-q0*q2),q0*q0-q1*q1-q2*q2+q3*q3);
  
  return;
}





/*
%
%	EP2Euler313(Q,E) translates the Euler parameter
%	vector Q into the corresponding (3-1-3) Euler angle
%	vector E.
%
*/
void	EP2Euler313(double *q, double *e)
{
  double t1, t2;
  
  t1 = atan2(q[3],q[2]);
  t2 = atan2(q[4],q[1]);

  e[1] = t1+t2;
  e[2] = 2*acos(sqrt(q[1]*q[1]+q[4]*q[4]));
  e[3] = t2-t1;
  
  return;
}




/*
%
%	EP2Euler321(Q,E) translates the Euler parameter vector
%	Q into the corresponding (3-2-1) Euler angle set.
%	
*/
void	EP2Euler321(double *q, double *e)
{
  double q0, q1, q2, q3;
  
  q0 = q[1];
  q1 = q[2];
  q2 = q[3];
  q3 = q[4];

  e[1] = atan2(2*(q1*q2+q0*q3),q0*q0+q1*q1-q2*q2-q3*q3);
  e[2] = asin(-2*(q1*q3-q0*q2));
  e[3]= atan2(2*(q2*q3+q0*q1),q0*q0-q1*q1-q2*q2+q3*q3);
  
  return;
}





/*
%
%	EP2Euler323(Q,E) translates the Euler parameter
%	vector Q into the corresponding (3-2-3) Euler angle
%	vector E.
%
*/
void	EP2Euler323(double *q, double *e)
{
  double t1, t2;

  t1 = atan2(q[2],q[3]);
  t2 = atan2(q[4],q[1]);

  e[1] = t2-t1;
  e[2] = 2*acos(sqrt(q[1]*q[1]+q[4]*q[4]));
  e[3] = t2+t1;
  
  return;
}




/*
%
%	EP2Gibbs(Q1,Q) translates the Euler parameter vector Q1
%	into the Gibbs vector Q.
%
*/
void	EP2Gibbs(double *q1, double *q)
{
  q[1] = q1[2]/q1[1];
  q[2] = q1[3]/q1[1];
  q[3] = q1[4]/q1[1];
  
  return;
}





/*
%
%	EP2MRP(Q1,Q) translates the Euler parameter vector Q1
%	into the MRP vector Q.
%
*/
void	EP2MRP(double *q1, double *q)
{
  q[1] = q1[2]/(1+q1[1]);
  q[2] = q1[3]/(1+q1[1]);
  q[3] = q1[4]/(1+q1[1]);
  
  return;
}




/*
%
%	EP2PRV(Q1,Q) translates the Euler parameter vector Q1
%	into the principal rotation vector Q.
%
*/
void	EP2PRV(double *q1, double *q)
{
  double p, sp;
  
  p = 2*acos(q1[1]);
  sp = sin(p/2);
  q[1] = q1[2]/sp*p;
  q[2] = q1[3]/sp*p;
  q[3] = q1[4]/sp*p;
  
  return;
}



/*
%   Euler1(X,M) 	Elementary rotation matrix
%	Returns the elementary rotation matrix about the
%	first body axis.
*/
void	Euler1(double x, double m[4][4])
{
  eye(m);
  m[2][2] = cos(x);
  m[2][3] = sin(x);
  m[3][2]= -m[2][3];
  m[3][3] = m[2][2];
  
  return;
}


/*
%   Euler2(X,M) 	Elementary rotation matrix
%	Returns the elementary rotation matrix about the
%	second body axis.
*/
void	Euler2(double x, double m[4][4])
{
  eye(m);
  m[1][1] = cos(x);
  m[1][3] = -sin(x);
  m[3][1]= -m[1][3];
  m[3][3] = m[1][1];
  
  return;
}



/*
%   Euler3(X,M) 	Elementary rotation matrix
%	Returns the elementary rotation matrix about the
%	third body axis.
*/
void	Euler3(double x, double m[4][4])
{
  eye(m);
  m[1][1] = cos(x);
  m[1][2] = sin(x);
  m[2][1]= -m[1][2];
  m[2][2] = m[1][1];
  
  return;
}




/*
%
%	Euler1212C(Q,C) returns the direction cosine 
%	matrix in terms of the 1-2-1 Euler angles.  
%	Input Q must be a 3x1 vector of Euler angles.
%	
*/
void	Euler1212C(double *q, double C[4][4])
{
  double st1, ct1, st2, ct2, st3, ct3;

  st1 = sin(q[1]);
  ct1 = cos(q[1]);
  st2 = sin(q[2]);
  ct2 = cos(q[2]);
  st3 = sin(q[3]);
  ct3 = cos(q[3]);

  C[1][1] = ct2;
  C[1][2] = st1*st2;
  C[1][3] = -ct1*st2;
  C[2][1] = st2*st3;
  C[2][2] = ct1*ct3-ct2*st1*st3;
  C[2][3] = ct3*st1+ct1*ct2*st3;
  C[3][1] = ct3*st2;
  C[3][2] = -ct2*ct3*st1-ct1*st3;
  C[3][3] = ct1*ct2*ct3-st1*st3;
  
  return;
}




/*
%
%	Euler1212EP(E,Q) translates the 121 Euler angle
%	vector E into the Euler parameter vector Q.
%
*/
void	Euler1212EP(double *e, double *q)
{
  double e1, e2, e3;
  
  e1 = e[1]/2;
  e2 = e[2]/2;
  e3 = e[3]/2;

  q[1] = cos(e2)*cos(e1+e3);
  q[2] = cos(e2)*sin(e1+e3);
  q[3] = sin(e2)*cos(e1-e3);
  q[4] = sin(e2)*sin(e1-e3);
  
  return;
}





/*
%
%	Euler1212Gibbs(E,Q) translates the (1-2-1) Euler
%	angle vector E into the Gibbs vector Q.
%
*/
void	Euler1212Gibbs(double *e, double *q)
{
  double ep[5];
  
  Euler1212EP(e,ep);
  EP2Gibbs(ep,q);
  
  return;
}





/*
%
%	Euler1212MRP(E,Q) translates the (1-2-1) Euler
%	angle vector E into the MRP vector Q.
%
*/
void	Euler1212MRP(double *e, double *q)
{
  double ep[5];
  
  Euler1212EP(e,ep);
  EP2MRP(ep,q);
  
  return;
}




/*
%
%	Euler1212PRV(E,Q) translates the (1-2-1) Euler
%	angle vector E into the principal rotation vector Q.
%
*/
void	Euler1212PRV(double *e, double *q)
{
  double ep[5];
  
  Euler1212EP(e,ep);
  EP2PRV(ep,q);
  
  return;
}




/*
%
%	Euler1232C(Q,C) returns the direction cosine 
%	matrix in terms of the 1-2-3 Euler angles.  
%	Input Q must be a 3x1 vector of Euler angles.
%
*/
void	Euler1232C(double *q, double C[4][4])
{
  double st1, st2, st3, ct1, ct2, ct3;
  
  st1 = sin(q[1]);
  ct1 = cos(q[1]);
  st2 = sin(q[2]);
  ct2 = cos(q[2]);
  st3 = sin(q[3]);
  ct3 = cos(q[3]);

  C[1][1] = ct2*ct3;
  C[1][2] = ct3*st1*st2+ct1*st3;
  C[1][3] = st1*st3-ct1*ct3*st2;
  C[2][1] = -ct2*st3;
  C[2][2] = ct1*ct3-st1*st2*st3;
  C[2][3] = ct3*st1+ct1*st2*st3;
  C[3][1] = st2;
  C[3][2] = -ct2*st1;
  C[3][3] = ct1*ct2;
  
  return;
}




/*
%
%	Euler1232EP(E,Q) translates the 123 Euler angle
%	vector E into the Euler parameter vector Q.
%
*/
void	Euler1232EP(double *e, double *q)
{
  double c1, c2, c3, s1, s2, s3;
  
  c1 = cos(e[1]/2);
  s1 = sin(e[1]/2);
  c2 = cos(e[2]/2);
  s2 = sin(e[2]/2);
  c3 = cos(e[3]/2);
  s3 = sin(e[3]/2);

  q[1] = c1*c2*c3-s1*s2*s3;
  q[2] = s1*c2*c3+c1*s2*s3;
  q[3] = c1*s2*c3-s1*c2*s3;
  q[4] = c1*c2*s3+s1*s2*c3;
  
  return;
}




/*
%
%	Euler1232Gibbs(E,Q) translates the (1-2-3) Euler
%	angle vector E into the Gibbs vector Q.
%
*/
void	Euler1232Gibbs(double *e, double *q)
{
  double ep[5];
  
  Euler1232EP(e,ep);
  EP2Gibbs(ep,q);
  
  return;
}





/*
%
%	Euler1232MRP(E,Q) translates the (1-2-3) Euler
%	angle vector E into the MRP vector Q.
%
*/
void	Euler1232MRP(double *e, double *q)
{
  double ep[5];
  
  Euler1232EP(e,ep);
  EP2MRP(ep,q);
  
  return;
}






/*
%
%	Euler1232PRV(E,Q) translates the (1-2-3) Euler
%	angle vector E into the principal rotation vector Q.
%
*/
void	Euler1232PRV(double *e, double *q)
{
  double ep[5];
  
  Euler1232EP(e,ep);
  EP2PRV(ep,q);
  
  return;
}




/*
%
%	Euler1312C(Q,C) returns the direction cosine 
%	matrix in terms of the 1-3-1 Euler angles.  
%	Input Q must be a 3x1 vector of Euler angles.
%	
*/
void	Euler1312C(double *q, double C[4][4])
{
  double st1, st2, st3, ct1, ct2, ct3;
  
  st1 = sin(q[1]);
  ct1 = cos(q[1]);
  st2 = sin(q[2]);
  ct2 = cos(q[2]);
  st3 = sin(q[3]);
  ct3 = cos(q[3]);

  C[1][1] = ct2;
  C[1][2] = ct1*st2;
  C[1][3] = st1*st2;
  C[2][1] = -ct3*st2;
  C[2][2] = ct1*ct2*ct3-st1*st3;
  C[2][3] = ct2*ct3*st1+ct1*st3;
  C[3][1] = st2*st3;
  C[3][2] = -ct3*st1-ct1*ct2*st3;
  C[3][3] = ct1*ct3-ct2*st1*st3;
  
  return;
}




/*
%
%	Euler1312EP(E,Q) translates the 131 Euler angle
%	vector E into the Euler parameter vector Q.
%
*/
void	Euler1312EP(double *e, double *q)
{
  double e1, e2, e3;
  
  e1 = e[1]/2;
  e2 = e[2]/2;
  e3 = e[3]/2;

  q[1] = cos(e2)*cos(e1+e3);
  q[2] = cos(e2)*sin(e1+e3);
  q[3] = sin(e2)*sin(-e1+e3);
  q[4] = sin(e2)*cos(-e1+e3);
  
  return;
}





/*
%
%	Euler1312Gibbs(E,Q) translates the (1-3-1) Euler
%	angle vector E into the Gibbs vector Q.
%
*/
void	Euler1312Gibbs(double *e, double *q)
{
  double ep[5];
  
  Euler1312EP(e,ep);
  EP2Gibbs(ep,q);
  
  return;
}




/*
%
%	Euler1312MRP(E,Q) translates the (1-3-1) Euler
%	angle vector E into the MRP vector Q.
%
*/
void	Euler1312MRP(double *e, double *q)
{
  double ep[5];
  
  Euler1312EP(e,ep);
  EP2MRP(ep,q);
  
  return;
}




/*
%
%	Euler1312PRV(E,Q) translates the (1-3-1) Euler
%	angle vector E into the principal rotation vector Q.
%
*/
void	Euler1312PRV(double *e, double *q)
{
  double ep[5];
  
  Euler1312EP(e,ep);
  EP2PRV(ep,q);
  
  return;
}



/*
%
%	Euler1322C(Q,C) returns the direction cosine 
%	matrix in terms of the 1-3-2 Euler angles.  
%	Input Q must be a 3x1 vector of Euler angles.
%
*/
void	Euler1322C(double *q, double C[4][4])
{
  double st1, st2, st3, ct1, ct2, ct3;
  
  st1 = sin(q[1]);
  ct1 = cos(q[1]);
  st2 = sin(q[2]);
  ct2 = cos(q[2]);
  st3 = sin(q[3]);
  ct3 = cos(q[3]);

  C[1][1] = ct2*ct3;
  C[1][2] = ct1*ct3*st2+st1*st3;
  C[1][3] = ct3*st1*st2-ct1*st3;
  C[2][1] = -st2;
  C[2][2] = ct1*ct2;
  C[2][3] = ct2*st1;
  C[3][1] = ct2*st3;
  C[3][2] = -ct3*st1+ct1*st2*st3;
  C[3][3] = ct1*ct3+st1*st2*st3;
  
  return;
}





/*
%
%	Euler1322EP(E,Q) translates the 132 Euler angle
%	vector E into the Euler parameter vector Q.
%
*/
void	Euler1322EP(double *e, double *q)
{
  double c1, c2, c3, s1, s2, s3;
  
  c1 = cos(e[1]/2);
  s1 = sin(e[1]/2);
  c2 = cos(e[2]/2);
  s2 = sin(e[2]/2);
  c3 = cos(e[3]/2);
  s3 = sin(e[3]/2);

  q[1] = c1*c2*c3+s1*s2*s3;
  q[2] = s1*c2*c3-c1*s2*s3;
  q[3] = c1*c2*s3-s1*s2*c3;
  q[4] = c1*s2*c3+s1*c2*s3;
  
  return;
}





/*
%
%	Euler1322Gibbs(E,Q) translates the (1-3-2) Euler
%	angle vector E into the Gibbs vector Q.
%
*/
void	Euler1322Gibbs(double *e, double *q)
{
  double ep[5];
  
  Euler1322EP(e,ep);
  EP2Gibbs(ep,q);
  
  return;
}





/*
%
%	Euler1322MRP(E,Q) translates the (1-3-2) Euler
%	angle vector E into the MRP vector Q.
%
*/
void	Euler1322MRP(double *e, double *q)
{
  double ep[5];
  
  Euler1322EP(e,ep);
  EP2MRP(ep,q);
  
  return;
}




/*
%
%	Euler1322PRV(E,Q) translates the (1-3-2) Euler
%	angle vector E into the principal rotation vector Q.
%
*/
void	Euler1322PRV(double *e, double *q)
{
  double ep[5];
  
  Euler1322EP(e,ep);
  EP2PRV(ep,q);
  
  return;
}




/*
%
%	Euler2122C(Q,C) returns the direction cosine 
%	matrix in terms of the 2-1-2 Euler angles.  
%	Input Q must be a 3x1 vector of Euler angles.
%	
*/
void	Euler2122C(double *q, double C[4][4])
{
  double st1, st2, st3, ct1, ct2, ct3;
  
  st1 = sin(q[1]);
  ct1 = cos(q[1]);
  st2 = sin(q[2]);
  ct2 = cos(q[2]);
  st3 = sin(q[3]);
  ct3 = cos(q[3]);

  C[1][1] = ct1*ct3-ct2*st1*st3;
  C[1][2] = st2*st3;
  C[1][3] = -ct3*st1-ct1*ct2*st3;
  C[2][1] = st1*st2;
  C[2][2] = ct2;
  C[2][3] = ct1*st2;
  C[3][1] = ct2*ct3*st1+ct1*st3;
  C[3][2] = -ct3*st2;
  C[3][3] = ct1*ct2*ct3-st1*st3;
  
  return;
}





/*
%
%	Euler2122EP(E,Q) translates the 212 Euler angle
%	vector E into the Euler parameter vector Q.
%
*/
void	Euler2122EP(double *e, double *q)
{
  double e1, e2, e3;
  
  e1 = e[1]/2;
  e2 = e[2]/2;
  e3 = e[3]/2;

  q[1] = cos(e2)*cos(e1+e3);
  q[2] = sin(e2)*cos(-e1+e3);
  q[3] = cos(e2)*sin(e1+e3);
  q[4] = sin(e2)*sin(-e1+e3);
  
  return;
}




/*
%
%	Euler2122Gibbs(E,Q) translates the (2-1-2) Euler
%	angle vector E into the Gibbs vector Q.
%
*/
void	Euler2122Gibbs(double *e, double *q)
{
  double ep[5];
  
  Euler2122EP(e,ep);
  EP2Gibbs(ep,q);
  
  return;
}





/*
%
%	Euler2122MRP(E,Q) translates the (2-1-2) Euler
%	angle vector E into the MRP vector Q.
%
*/
void	Euler2122MRP(double *e, double *q)
{
  double ep[5];
  
  Euler2122EP(e,ep);
  EP2MRP(ep,q);
  
  return;
}




/*
%
%	Euler2122PRV(E,Q) translates the (2-1-2) Euler
%	angle vector E into the principal rotation vector Q.
%
*/
void	Euler2122PRV(double *e, double *q)
{
  double ep[5];
  
  Euler2122EP(e,ep);
  EP2PRV(ep,q);
  
  return;
}




/*
%
%	Euler2132C(Q,C) returns the direction cosine 
%	matrix in terms of the 2-1-3 Euler angles.  
%	Input Q must be a 3x1 vector of Euler angles.
%	
*/
void	Euler2132C(double *q, double C[4][4])
{
  double st1, st2, st3, ct1, ct2, ct3;
  
  st1 = sin(q[1]);
  ct1 = cos(q[1]);
  st2 = sin(q[2]);
  ct2 = cos(q[2]);
  st3 = sin(q[3]);
  ct3 = cos(q[3]);

  C[1][1] = ct1*ct3+st1*st2*st3;
  C[1][2] = ct2*st3;
  C[1][3] = -ct3*st1+ct1*st2*st3;
  C[2][1] = ct3*st1*st2-ct1*st3;
  C[2][2] = ct2*ct3;
  C[2][3] = ct1*ct3*st2 + st1*st3;
  C[3][1] = ct2*st1;
  C[3][2] = -st2;
  C[3][3] = ct1*ct2;
  
  return;
}




/*
%
%	Euler2132EP(E,Q) translates the 213 Euler angle
%	vector E into the Euler parameter vector Q.
%
*/
void	Euler2132EP(double *e, double *q)
{
  double c1, c2, c3, s1, s2, s3;
  
  c1 = cos(e[1]/2);
  s1 = sin(e[1]/2);
  c2 = cos(e[2]/2);
  s2 = sin(e[2]/2);
  c3 = cos(e[3]/2);
  s3 = sin(e[3]/2);

  q[1] = c1*c2*c3+s1*s2*s3;
  q[2] = c1*s2*c3+s1*c2*s3;
  q[3] = s1*c2*c3-c1*s2*s3;
  q[4] = c1*c2*s3-s1*s2*c3;
  
  return;
}





/*
%
%	Euler2132Gibbs(E,Q) translates the (2-1-3) Euler
%	angle vector E into the Gibbs vector Q.
%
*/
void	Euler2132Gibbs(double *e, double *q)
{
  double ep[5];
  
  Euler2132EP(e,ep);
  EP2Gibbs(ep,q);
  
  return;
}





/*
%
%	Euler2132MRP(E,Q) translates the (2-1-3) Euler
%	angle vector E into the MRP vector Q.
%
*/
void	Euler2132MRP(double *e, double *q)
{
  double ep[5];
  
  Euler2132EP(e,ep);
  EP2MRP(ep,q);
  
  return;
}




/*
%
%	Euler2132PRV(E,Q) translates the (2-1-3) Euler
%	angle vector E into the principal rotation vector Q.
%
*/
void	Euler2132PRV(double *e, double *q)
{
  double ep[5];
  
  Euler2132EP(e,ep);
  EP2PRV(ep,q);
  
  return;
}





/*
%
%	Euler2312C(Q,C) returns the direction cosine 
%	matrix in terms of the 2-3-1 Euler angles.  
%	Input Q must be a 3x1 vector of Euler angles.
%	
*/
void	Euler2312C(double *q, double C[4][4])
{
  double st1, st2, st3, ct1, ct2, ct3;
  
  st1 = sin(q[1]);
  ct1 = cos(q[1]);
  st2 = sin(q[2]);
  ct2 = cos(q[2]);
  st3 = sin(q[3]);
  ct3 = cos(q[3]);

  C[1][1] = ct1*ct2;
  C[1][2] = st2;
  C[1][3] = -ct2*st1;
  C[2][1] = -ct1*ct3*st2+st1*st3;
  C[2][2] = ct2*ct3;
  C[2][3] = ct3*st1*st2+ct1*st3;
  C[3][1]= ct3*st1+ct1*st2*st3;
  C[3][2] = -ct2*st3;
  C[3][3] = ct1*ct3-st1*st2*st3;
  
  return;
}





/*
%
%	Euler2312EP(E,Q) translates the 231 Euler angle
%	vector E into the Euler parameter vector Q.
%
*/
void	Euler2312EP(double *e, double *q)
{
  double c1, c2, c3, s1, s2, s3;
  
  c1 = cos(e[1]/2);
  s1 = sin(e[1]/2);
  c2 = cos(e[2]/2);
  s2 = sin(e[2]/2);
  c3 = cos(e[3]/2);
  s3 = sin(e[3]/2);

  q[1] = c1*c2*c3-s1*s2*s3;
  q[2] = c1*c2*s3+s1*s2*c3;
  q[3] = s1*c2*c3+c1*s2*s3;
  q[4] = c1*s2*c3-s1*c2*s3;
  
  return;
}






/*
%
%	Euler2312Gibbs(E,Q) translates the (2-3-1) Euler
%	angle vector E into the Gibbs vector Q.
%
*/
void	Euler2312Gibbs(double *e, double *q)
{
  double ep[5];
  
  Euler2312EP(e,ep);
  EP2Gibbs(ep,q);
  
  return;
}




/*
%
%	Euler2312MRP(E,Q) translates the (2-3-1) Euler
%	angle vector E into the MRP vector Q.
%
*/
void	Euler2312MRP(double *e, double *q)
{
  double ep[5];
  
  Euler2312EP(e,ep);
  EP2MRP(ep,q);
  
  return;
}





/*
%
%	Euler2312PRV(E,Q) translates the (2-3-1) Euler
%	angle vector E into the principal rotation vector Q.
%
*/
void	Euler2312PRV(double *e, double *q)
{
  double ep[5];
  
  Euler2312EP(e,ep);
  EP2PRV(ep,q);
  
  return;
}





/*
%
%	Euler2322C(Q) returns the direction cosine 
%	matrix in terms of the 2-3-2 Euler angles.  
%	Input Q must be a 3x1 vector of Euler angles.
%	
*/
void	Euler2322C(double *q, double C[4][4])
{
  double st1, st2, st3, ct1, ct2, ct3;
  
  st1 = sin(q[1]);
  ct1 = cos(q[1]);
  st2 = sin(q[2]);
  ct2 = cos(q[2]);
  st3 = sin(q[3]);
  ct3 = cos(q[3]);

  C[1][1] = ct1*ct2*ct3-st1*st3;
  C[1][2] = ct3*st2;
  C[1][3] = -ct2*ct3*st1-ct1*st3;
  C[2][1] = -ct1*st2;
  C[2][2] = ct2;
  C[2][3] = st1*st2;
  C[3][1] = ct3*st1+ct1*ct2*st3;
  C[3][2] = st2*st3;
  C[3][3] = ct1*ct3-ct2*st1*st3;
  
  return;
}




/*
%
%	Euler2322EP(E,Q) translates the 232 Euler angle
%	vector E into the Euler parameter vector Q.
%
*/
void	Euler2322EP(double *e, double *q)
{
  double e1, e2, e3;
  
  e1 = e[1]/2;
  e2 = e[2]/2;
  e3 = e[3]/2;

  q[1] = cos(e2)*cos(e1+e3);
  q[2] = sin(e2)*sin(e1-e3);
  q[3] = cos(e2)*sin(e1+e3);
  q[4] = sin(e2)*cos(e1-e3);
  
  return;
}




/*
%
%	Euler2322Gibbs(E) translates the (2-3-2) Euler
%	angle vector E into the Gibbs vector Q.
%
*/
void	Euler2322Gibbs(double *e, double *q)
{
  double ep[5];
  
  Euler2322EP(e,ep);
  EP2Gibbs(ep,q);
  
  return;
}





/*
%
%	Euler2322MRP(E,Q) translates the (2-3-2) Euler
%	angle vector E into the MRP vector Q.
%
*/
void	Euler2322MRP(double *e, double *q)
{
  double ep[5];
  
  Euler2322EP(e,ep);
  EP2MRP(ep,q);
  
  return;
}




/*
%
%	Euler2322PRV(E,Q) translates the (2-3-2) Euler
%	angle vector E into the principal rotation vector Q.
%
*/
void	Euler2322PRV(double *e, double *q)
{
  double ep[5];
  
  Euler2322EP(e,ep);
  EP2PRV(ep,q);
  
  return;
}




/*
%
%	Euler3122C(Q,C) returns the direction cosine 
%	matrix in terms of the 1-2-3 Euler angles.  
%	Input Q must be a 3x1 vector of Euler angles.
%
*/
void	Euler3122C(double *q, double C[4][4])
{
  double st1, st2, st3, ct1, ct2, ct3;
  
  st1 = sin(q[1]);
  ct1 = cos(q[1]);
  st2 = sin(q[2]);
  ct2 = cos(q[2]);
  st3 = sin(q[3]);
  ct3 = cos(q[3]);

  C[1][1] = ct1*ct3-st1*st2*st3;
  C[1][2] = ct3*st1+ct1*st2*st3;
  C[1][3] = -ct2*st3;
  C[2][1] = -ct2*st1;
  C[2][2] = ct1*ct2;
  C[2][3] = st2;
  C[3][1] = ct3*st1*st2+ct1*st3;
  C[3][2] = st1*st3-ct1*ct3*st2;
  C[3][3] = ct2*ct3;
  
  return;
}




/*
%
%	Euler3122EP(E,Q) translates the 312 Euler angle
%	vector E into the Euler parameter vector Q.
%
*/
void	Euler3122EP(double *e, double *q)
{
  double c1, c2, c3, s1, s2, s3;
  
  c1 = cos(e[1]/2);
  s1 = sin(e[1]/2);
  c2 = cos(e[2]/2);
  s2 = sin(e[2]/2);
  c3 = cos(e[3]/2);
  s3 = sin(e[3]/2);

  q[1] = c1*c2*c3-s1*s2*s3;
  q[2] = c1*s2*c3-s1*c2*s3;
  q[3] = c1*c2*s3+s1*s2*c3;
  q[4] = s1*c2*c3+c1*s2*s3;
  
  return;
}







/*
%
%	Euler3122Gibbs(E,Q) translates the (3-1-2) Euler
%	angle vector E into the Gibbs vector Q.
%
*/
void	Euler3122Gibbs(double *e, double *q)
{
  double ep[5];
  
  Euler3122EP(e,ep);
  EP2Gibbs(ep,q);
  
  return;
}




/*
%
%	Euler3122MRP(E,Q) translates the (3-1-2) Euler
%	angle vector E into the MRP vector Q.
%
*/
void	Euler3122MRP(double *e, double *q)
{
  double ep[5];
  
  Euler3122EP(e,ep);
  EP2MRP(ep,q);
  
  return;
}




/*
%
%	Euler3122PRV(E,Q) translates the (3-1-2) Euler
%	angle vector E into the principal rotation vector Q.
%
*/
void	Euler3122PRV(double *e, double *q)
{
  double ep[5];
  
  Euler3122EP(e,ep);
  EP2PRV(ep,q);
  
  return;
}




/*
%
%	Euler3132C(Q,C) returns the direction cosine 
%	matrix in terms of the 3-1-3 Euler angles.  
%	Input Q must be a 3x1 vector of Euler angles.
%
*/
void	Euler3132C(double *q, double C[4][4])
{
  double st1, st2, st3, ct1, ct2, ct3;
  
  st1 = sin(q[1]);
  ct1 = cos(q[1]);
  st2 = sin(q[2]);
  ct2 = cos(q[2]);
  st3 = sin(q[3]);
  ct3 = cos(q[3]);

  C[1][1] = ct3*ct1-st3*ct2*st1;
  C[1][2] = ct3*st1+st3*ct2*ct1;
  C[1][3] = st3*st2;
  C[2][1] = -st3*ct1-ct3*ct2*st1;
  C[2][2] = -st3*st1+ct3*ct2*ct1;
  C[2][3] = ct3*st2;
  C[3][1] = st2*st1;
  C[3][2] = -st2*ct1;
  C[3][3] = ct2;
  
  return;
}



/*
%
%	Euler3132EP(E,Q) translates the 313 Euler angle
%	vector E into the Euler parameter vector Q.
%
*/
void	Euler3132EP(double *e, double *q)
{
  double e1, e2, e3;
  
  e1 = e[1]/2;
  e2 = e[2]/2;
  e3 = e[3]/2;

  q[1] = cos(e2)*cos(e1+e3);
  q[2] = sin(e2)*cos(e1-e3);
  q[3] = sin(e2)*sin(e1-e3);
  q[4] = cos(e2)*sin(e1+e3);
  
  return;
}






/*
%
%	Euler3132Gibbs(E,Q) translates the (3-1-3) Euler
%	angle vector E into the Gibbs vector Q.
%
*/
void	Euler3132Gibbs(double *e, double *q)
{
  double ep[5];
  
  Euler3132EP(e,ep);
  EP2Gibbs(ep,q);
  
  return;
}





/*
%
%	Euler3132MRP(E,Q) translates the (3-1-3) Euler
%	angle vector E into the MRP vector Q.
%
*/
void	Euler3132MRP(double *e, double *q)
{
  double ep[5];
  
  Euler3132EP(e,ep);
  EP2MRP(ep,q);
  
  return;
}




/*
%
%	Euler3132PRV(E,Q) translates the (3-1-3) Euler
%	angle vector E into the principal rotation vector Q.
%
*/
void	Euler3132PRV(double *e, double *q)
{
  double ep[5];
  
  Euler3132EP(e,ep);
  EP2PRV(ep,q);
  
  return;
}




/*
%
%	Euler3212C(Q,C) returns the direction cosine 
%	matrix in terms of the 3-2-1 Euler angles.  
%	Input Q must be a 3x1 vector of Euler angles.
%
*/
void	Euler3212C(double *q, double C[4][4])
{
  double st1, st2, st3, ct1, ct2, ct3;
  
  st1 = sin(q[1]);
  ct1 = cos(q[1]);
  st2 = sin(q[2]);
  ct2 = cos(q[2]);
  st3 = sin(q[3]);
  ct3 = cos(q[3]);

  C[1][1] = ct2*ct1;
  C[1][2] = ct2*st1;
  C[1][3] = -st2;
  C[2][1] = st3*st2*ct1-ct3*st1;
  C[2][2] = st3*st2*st1+ct3*ct1;
  C[2][3] = st3*ct2;
  C[3][1] = ct3*st2*ct1+st3*st1;
  C[3][2] = ct3*st2*st1-st3*ct1;
  C[3][3] = ct3*ct2;
  
  return;
}





/*
%
%	Euler3212EPE,Q) translates the 321 Euler angle
%	vector E into the Euler parameter vector Q.
%
*/
void	Euler3212EP(double *e, double *q)
{
  double c1, c2, c3, s1, s2, s3;
  
  c1 = cos(e[1]/2);
  s1 = sin(e[1]/2);
  c2 = cos(e[2]/2);
  s2 = sin(e[2]/2);
  c3 = cos(e[3]/2);
  s3 = sin(e[3]/2);

  q[1] = c1*c2*c3+s1*s2*s3;
  q[2] = c1*c2*s3-s1*s2*c3;
  q[3] = c1*s2*c3+s1*c2*s3;
  q[4] = s1*c2*c3-c1*s2*s3;
  
  return;
}




/*
%
%	Euler3212Gibbs(E,Q) translates the (3-2-1) Euler
%	angle vector E into the Gibbs vector Q.
%
*/
void	Euler3212Gibbs(double *e, double *q)
{
  double ep[5];
  
  Euler3212EP(e,ep);
  EP2Gibbs(ep,q);
  
  return;
}




/*
%
%	Euler3212MRP(E,Q) translates the (3-2-1) Euler
%	angle vector E into the MRP vector Q.
%
*/
void	Euler3212MRP(double *e, double *q)
{
  double ep[5];
  
  Euler3212EP(e,ep);
  EP2MRP(ep,q);
  
  return;
}




/*
%
%	Euler3212PRV(E,Q) translates the (3-2-1) Euler
%	angle vector E into the principal rotation vector Q.
%
*/
void	Euler3212PRV(double *e, double *q)
{
  double ep[5];
  
  Euler3212EP(e,ep);
  EP2PRV(ep,q);
  
  return;
}




/*
%
%	Euler3232C(Q,C) returns the direction cosine 
%	matrix in terms of the 3-2-3 Euler angles.  
%	Input Q must be a 3x1 vector of Euler angles.
%	
*/
void	Euler3232C(double *q, double C[4][4])
{
  double st1, st2, st3, ct1, ct2, ct3;
  
  st1 = sin(q[1]);
  ct1 = cos(q[1]);
  st2 = sin(q[2]);
  ct2 = cos(q[2]);
  st3 = sin(q[3]);
  ct3 = cos(q[3]);

  C[1][1] = ct1*ct2*ct3-st1*st3;
  C[1][2] = ct2*ct3*st1+ct1*st3;
  C[1][3] = -ct3*st2;
  C[2][1] = -ct3*st1-ct1*ct2*st3;
  C[2][2] = ct1*ct3-ct2*st1*st3;
  C[2][3] = st2*st3;
  C[3][1] = ct1*st2;
  C[3][2] = st1*st2;
  C[3][3] = ct2;
  
  return;
}





/*
%
%	Euler3232EP(E,Q) translates the 323 Euler angle
%	vector E into the Euler parameter vector Q.
%
*/
void	Euler3232EP(double *e, double *q)
{
  double e1, e2, e3;
  
  e1 = e[1]/2;
  e2 = e[2]/2;
  e3 = e[3]/2;

  q[1] = cos(e2)*cos(e1+e3);
  q[2] = sin(e2)*sin(-e1+e3);
  q[3] = sin(e2)*cos(-e1+e3);
  q[4] = cos(e2)*sin(e1+e3);
  
  return;
}







/*
%
%	Euler3232Gibbs(E,Q) translates the (3-2-3) Euler
%	angle vector E into the Gibbs vector Q.
%
*/
void	Euler3232Gibbs(double *e, double *q)
{
  double ep[5];
  
  Euler3232EP(e,ep);
  EP2Gibbs(ep,q);
  
  return;
}





/*
%
%	Euler3232MRP(E,Q) translates the (3-2-3) Euler
%	angle vector E into the MRP vector Q.
%
*/
void	Euler3232MRP(double *e, double *q)
{
  double ep[5];
  
  Euler3232EP(e,ep);
  EP2MRP(ep,q);
  
  return;
}





/*
%
%	Euler3232PRV(E,Q) translates the (3-2-3) Euler
%	angle vector Q1 into the principal rotation vector Q.
%
*/
void	Euler3232PRV(double *e, double *q)
{
  double ep[5];
  
  Euler3232EP(e,ep);
  EP2PRV(ep,q);
  
  return;
}




/*
%
%	Gibbs2C(Q,C) returns the direction cosine 
%	matrix in terms of the 3x1 Gibbs vector Q.  
%	
*/
void	Gibbs2C(double *q, double C[4][4])
{
  double q1, q2, q3, d1;
  
  q1 = q[1];
  q2 = q[2];
  q3 = q[3];

  d1 = dot(q,q);
  C[1][1] = 1+2*q1*q1-d1;
  C[1][2] = 2*(q1*q2+q3);
  C[1][3] = 2*(q1*q3-q2);
  C[2][1] = 2*(q2*q1-q3);
  C[2][2] = 1+2*q2*q2-d1;
  C[2][3] = 2*(q2*q3+q1);
  C[3][1] = 2*(q3*q1+q2);
  C[3][2] = 2*(q3*q2-q1);
  C[3][3] = 1+2*q3*q3-d1;
  Mmult( 1./(1+d1),C, C);
  
  return;
}





/*
%
%	Gibbs2EP(Q1,Q) translates the Gibbs vector Q1
%	into the Euler parameter vector Q.
%
*/
void	Gibbs2EP(double *q1, double *q)
{

  q[1] = 1/sqrt(1+dot(q1,q1));
  q[2] = q1[1]*q[1];
  q[3] = q1[2]*q[1];
  q[4] = q1[3]*q[1];
  
  return;
}




/*
%
%	 Gibbs2Euler121(Q,E) translates the Gibbs
%	 vector Q into the (1-2-1) Euler angle vector E.
%
*/
void	Gibbs2Euler121(double *q, double *e)
{
  double ep[5];
  
  Gibbs2EP(q,ep);
  EP2Euler121(ep,e);
  
  return;
}




/*
%
%	 Gibbs2Euler123(Q,E) translates the Gibbs
%	 vector Q into the (1-2-3) Euler angle vector E.
%
*/
void	Gibbs2Euler123(double *q, double *e)
{
  double ep[5];
  
  Gibbs2EP(q,ep);
  EP2Euler123(ep,e);
  
  return;
}




/*
%
%	 Gibbs2Euler131(Q,E) translates the Gibbs
%	 vector Q into the (1-3-1) Euler angle vector E.
%
*/
void	Gibbs2Euler131(double *q, double *e)
{
  double ep[5];
  
  Gibbs2EP(q,ep);
  EP2Euler131(ep,e);
  
  return;
}





/*
%
%	 Gibbs2Euler132(Q,E) translates the Gibbs
%	 vector Q into the (1-3-2) Euler angle vector E.
%
*/
void	Gibbs2Euler132(double *q, double *e)
{
  double ep[5];
  
  Gibbs2EP(q,ep);
  EP2Euler132(ep,e);
  
  return;
}





/*
%
%	 Gibbs2Euler212(Q,E) translates the Gibbs
%	 vector Q into the (2-1-2) Euler angle vector E.
%
*/
void	Gibbs2Euler212(double *q, double *e)
{
  double ep[5];
  
  Gibbs2EP(q,ep);
  EP2Euler212(ep,e);
  
  return;
}





/*
%
%	 Gibbs2Euler213(Q,E) translates the Gibbs
%	 vector Q into the (2-1-3) Euler angle vector E.
%
*/
void	Gibbs2Euler213(double *q, double *e)
{
  double ep[5];
  
  Gibbs2EP(q,ep);
  EP2Euler213(ep,e);
  
  return;
}





/*
%
%	 Gibbs2Euler231(Q,E) translates the Gibbs
%	 vector Q into the (2-3-1) Euler angle vector E.
%
*/
void	Gibbs2Euler231(double *q, double *e)
{
  double ep[5];
  
  Gibbs2EP(q,ep);
  EP2Euler231(ep,e);
  
  return;
}





/*
%
%	 Gibbs2Euler232(Q,E) translates the Gibbs
%	 vector Q into the (2-3-2) Euler angle vector E.
%
*/
void	Gibbs2Euler232(double *q, double *e)
{
  double ep[5];
  
  Gibbs2EP(q,ep);
  EP2Euler232(ep,e);
  
  return;
}





/*
%
%	 Gibbs2Euler312(Q,E) translates the Gibbs
%	 vector Q into the (3-1-2) Euler angle vector E.
%
*/
void	Gibbs2Euler312(double *q, double *e)
{
  double ep[5];
  
  Gibbs2EP(q,ep);
  EP2Euler312(ep,e);
  
  return;
}





/*
%
%	 Gibbs2Euler313(Q,E) translates the Gibbs
%	 vector Q into the (3-1-3) Euler angle vector E.
%
*/
void	Gibbs2Euler313(double *q, double *e)
{
  double ep[5];
  
  Gibbs2EP(q,ep);
  EP2Euler313(ep,e);
  
  return;
}





/*
%
%	 Gibbs2Euler321(Q,E) translates the Gibbs
%	 vector Q into the (3-2-1) Euler angle vector E.
%
*/
void	Gibbs2Euler321(double *q, double *e)
{
  double ep[5];
  
  Gibbs2EP(q,ep);
  EP2Euler321(ep,e);
  
  return;
}





/*
%
%	 Gibbs2Euler323(Q,E) translates the Gibbs
%	 vector Q into the (3-2-3) Euler angle vector E.
%
*/
void	Gibbs2Euler323(double *q, double *e)
{
  double ep[5];
  
  Gibbs2EP(q,ep);
  EP2Euler323(ep,e);
  
  return;
}






/*
%
%	Gibbs2MRP(Q1,Q) translates the Gibbs vector Q1
%	into the MRP vector Q.
%
*/
void	Gibbs2MRP(double *q1, double *q)
{
  
  mult(1.0/(1+sqrt(1+dot(q1,q1))),q1,q);
  
  return;
}




/*
%
%	Gibbs2PRV(Q1,Q) translates the Gibbs vector Q1
%	into the principal rotation vector Q.
%
*/
void	Gibbs2PRV(double *q1, double *q)
{
  double tp, p;
  
  tp = sqrt(dot(q1,q1));
  p = 2*atan(tp);
  q[1] = q1[1]/tp*p;
  q[2] = q1[2]/tp*p;
  q[3] = q1[3]/tp*p;
  
  return;
}




/*
%
%	MRP2C(Q,C) returns the direction cosine 
%	matrix in terms of the 3x1 MRP vector Q.  
%	
*/
void	MRP2C(double *q, double C[4][4])
{
  double q1, q2, q3, S, d1, d;
  
  q1 = q[1];
  q2 = q[2];
  q3 = q[3];

  d1 = dot(q,q);
  S = 1-d1;
  d = (1+d1)*(1+d1);
  C[1][1] = 4*(2*q1*q1-d1)+S*S;
  C[1][2] = 8*q1*q2+4*q3*S;
  C[1][3] = 8*q1*q3-4*q2*S;
  C[2][1] = 8*q2*q1-4*q3*S;
  C[2][2] = 4*(2*q2*q2-d1)+S*S;
  C[2][3] = 8*q2*q3+4*q1*S;
  C[3][1] = 8*q3*q1+4*q2*S;
  C[3][2] = 8*q3*q2-4*q1*S;
  C[3][3] = 4*(2*q3*q3-d1)+S*S;
  Mmult(1./d, C, C);
  
  return;
}




/*
%
%	MRP2EP(Q1,Q) translates the MRP vector Q1
%	into the Euler parameter vector Q.
%
*/
void	MRP2EP(double *q1, double *q)
{
  double ps;
  
  ps = 1+dot(q1,q1);
  q[1] = (1-dot(q1,q1))/ps;
  q[2] = 2*q1[1]/ps;
  q[3] = 2*q1[2]/ps;
  q[4] = 2*q1[3]/ps;
  
  return;
}




/*
%
%	 MRP2Euler121(Q,E) translates the MRP
%	 vector Q into the (1-2-1) Euler angle vector E.
%
*/
void	MRP2Euler121(double *q, double *e)
{
  double ep[5];
  
  MRP2EP(q,ep);
  EP2Euler121(ep,e);
  
  return;
}






/*
%
%	 MRP2Euler123(Q,E) translates the MRP
%	 vector Q into the (1-2-3) Euler angle vector E.
%
*/
void	MRP2Euler123(double *q, double *e)
{
  double ep[5];
  
  MRP2EP(q,ep);
  EP2Euler123(ep,e);
  
  return;
}




/*
%
%	 MRP2Euler131(Q,E) translates the MRP
%	 vector Q into the (1-3-1) Euler angle vector E.
%
*/
void	MRP2Euler131(double *q, double *e)
{
  double ep[5];
  
  MRP2EP(q,ep);
  EP2Euler131(ep,e);
  
  return;
}




/*
%
%	 MRP2Euler132(Q,E) translates the MRP
%	 vector Q into the (1-3-2) Euler angle vector E.
%
*/
void	MRP2Euler132(double *q, double *e)
{
  double ep[5];
  
  MRP2EP(q,ep);
  EP2Euler132(ep,e);
  
  return;
}




/*
%
%	E = MRP2Euler212(Q) translates the MRP
%	 vector Q into the (2-1-2) Euler angle vector E.
%
*/
void	MRP2Euler212(double *q, double *e)
{
  double ep[5];
  
  MRP2EP(q,ep);
  EP2Euler212(ep,e);
  
  return;
}




/*
%
%	 MRP2Euler213(Q,E) translates the MRP
%	 vector Q into the (2-1-3) Euler angle vector E.
%
*/
void	MRP2Euler213(double *q, double *e)
{
  double ep[5];
  
  MRP2EP(q,ep);
  EP2Euler213(ep,e);
  
  return;
}





/*
%
%	 MRP2Euler231(Q,E) translates the MRP
%	 vector Q into the (2-3-1) Euler angle vector E.
%
*/
void	MRP2Euler231(double *q, double *e)
{
  double ep[5];
  
  MRP2EP(q,ep);
  EP2Euler231(ep,e);
  
  return;
}




/*
%
%	 MRP2Euler232(Q,E) translates the MRP
%	 vector Q into the (2-3-2) Euler angle vector E.
%
*/
void	MRP2Euler232(double *q, double *e)
{
  double ep[5];
  
  MRP2EP(q,ep);
  EP2Euler232(ep,e);
  
  return;
}




/*
%
%	 MRP2Euler312(Q,E) translates the MRP
%	 vector Q into the (3-1-2) Euler angle vector E.
%
*/
void	MRP2Euler312(double *q, double *e)
{
  double ep[5];
  
  MRP2EP(q,ep);
  EP2Euler312(ep,e);
  
  return;
}




/*
%
%	 MRP2Euler313(Q,E) translates the MRP
%	 vector Q into the (3-1-3) Euler angle vector E.
%
*/
void	MRP2Euler313(double *q, double *e)
{
  double ep[5];
  
  MRP2EP(q,ep);
  EP2Euler313(ep,e);
  
  return;
}





/*
%
%	 MRP2Euler321(Q,E) translates the MRP
%	 vector Q into the (3-2-1) Euler angle vector E.
%
*/
void	MRP2Euler321(double *q, double *e)
{
  double ep[5];
  
  MRP2EP(q,ep);
  EP2Euler321(ep,e);
  
  return;
}





/*
%
%	 MRP2Euler323(Q,E) translates the MRP
%	 vector Q into the (3-2-3) Euler angle vector E.
%
*/
void	MRP2Euler323(double *q, double *e)
{
  double ep[5];
  
  MRP2EP(q,ep);
  EP2Euler323(ep,e);
  
  return;
}



/*
%
%	MRP2Gibbs(Q1,Q) translates the MRP vector Q1
%	into the Gibbs vector Q.
%
*/
void	MRP2Gibbs(double *q1, double *q)
{
  mult(2./(1.-dot(q1,q1)),q1, q);
  
  return;
}




/*
%
%	MRP2PRV(Q1,Q) translates the MRP vector Q1
%	into the principal rotation vector Q.
%
*/
void	MRP2PRV(double *q1, double *q)
{
  double tp, p;
  
  tp = sqrt(dot(q1,q1));
  p = 4*atan(tp);
  q[1] = q1[1]/tp*p;
  q[2] = q1[2]/tp*p;
  q[3] = q1[3]/tp*p;
  
  return;
}






/*
%
%	MRPswitch(Q,s2,s) checks to see if norm(Q) is larger than s2.
%	If yes, then the MRP vector Q is mapped to its shadow set.
%
*/
void	MRPswitch(double *q, double s2, double *s)
{
  double q2;
  
  q2 = dot(q,q);
  if (q2>s2*s2) {
	  mult(-1./q2,q,s);
  } else {
	  equal(q,s);
  }
  
  return;
}





/*
%
% 	Makes sure that the angle x lies within +/- Pi.
%
*/
double	Picheck(double x)
{
  double q;
  
  q = x;

  if (x >  M_PI) q = x-2*M_PI;
  
  if (x < -M_PI) q = x+2*M_PI;
  
  
  return q;
}





/*
%
%	PRV2C(Q,C) returns the direction cosine 
%	matrix in terms of the 3x1 principal rotation vector
%	Q.  
%
*/
void	PRV2C(double *q, double C[4][4])
{
  double q0, q1, q2, q3, cp, sp, d1;
  
  q0 = sqrt(dot(q,q));
  q1 = q[1]/q0;
  q2 = q[2]/q0;
  q3 = q[3]/q0;

  cp= cos(q0);
  sp= sin(q0);
  d1 = 1-cp;
  C[1][1] = q1*q1*d1+cp;
  C[1][2] = q1*q2*d1+q3*sp;
  C[1][3] = q1*q3*d1-q2*sp;
  C[2][1] = q2*q1*d1-q3*sp;
  C[2][2] = q2*q2*d1+cp;
  C[2][3] = q2*q3*d1+q1*sp;
  C[3][1] = q3*q1*d1+q2*sp;
  C[3][2] = q3*q2*d1-q1*sp;
  C[3][3] = q3*q3*d1+cp;
  
  return;
}




/*
%	
%	PRV2elem(R,Q) translates a prinicpal rotation vector R
%	into the corresponding principal rotation element set Q.
%
*/
void	PRV2elem(double *r, double *q)
{
  
  q[1] = sqrt(dot(r,r));
  q[2] = r[1]/q[1];
  q[3] = r[2]/q[1];
  q[4] = r[3]/q[1];
  
  return;
}






/*
%
%	PRV2EP(Q0,Q) translates the principal rotation vector Q1
%	into the Euler parameter vector Q.
%
*/
void	PRV2EP(double *q0, double *q)
{
  double q1[5], sp;
  
  PRV2elem(q0,q1);
  sp = sin(q1[1]/2);
  q[1] = cos(q1[1]/2);
  q[2] = q1[2]*sp;
  q[3] = q1[3]*sp;
  q[4] = q1[4]*sp;
  
  return;
}





/*
%
%	PRV2Euler121(Q,E) translates the principal rotation
%	vector Q into the (1-2-1) Euler angle vector E.
%
*/
void	PRV2Euler121(double *q, double *e)
{
  double ep[5];
  
  PRV2EP(q,ep);
  EP2Euler121(ep,e);
  
  return;
}





/*
%
%	PRV2Euler123(Q,E) translates the principal rotation
%	vector Q into the (1-2-3) Euler angle vector E.
%
*/
void	PRV2Euler123(double *q, double *e)
{
  double ep[5];
  
  PRV2EP(q,ep);
  EP2Euler123(ep,e);
  
  return;
}





/*
%
%	PRV2Euler131(Q,E) translates the principal rotation
%	vector Q into the (1-3-1) Euler angle vector E.
%
*/
void	PRV2Euler131(double *q, double *e)
{
  double ep[5];
  
  PRV2EP(q,ep);
  EP2Euler131(ep,e);
  
  return;
}





/*
%
%	PRV2Euler132(Q,E) translates the principal rotation
%	vector Q into the (1-3-2) Euler angle vector E.
%
*/
void	PRV2Euler132(double *q, double *e)
{
  double ep[5];
  
  PRV2EP(q,ep);
  EP2Euler132(ep,e);
  
  return;
}






/*
%
%	PRV2Euler212(Q,E) translates the principal rotation
%	vector Q into the (2-1-2) Euler angle vector E.
%
*/
void	PRV2Euler212(double *q, double *e)
{
  double ep[5];
  
  PRV2EP(q,ep);
  EP2Euler212(ep,e);
  
  return;
}





/*
%
%	PRV2Euler213(Q,E) translates the principal rotation
%	vector Q into the (2-1-3) Euler angle vector E.
%
*/
void	PRV2Euler213(double *q, double *e)
{
  double ep[5];
  
  PRV2EP(q,ep);
  EP2Euler213(ep,e);
  
  return;
}





/*
%
%	PRV2Euler231(Q) translates the principal rotation
%	vector Q into the (2-3-1) Euler angle vector E.
%
*/
void	PRV2Euler231(double *q, double *e)
{
  double ep[5];
  
  PRV2EP(q,ep);
  EP2Euler231(ep,e);
  
  return;
}





/*
%
%	PRV2Euler232(Q,E) translates the principal rotation
%	vector Q into the (2-3-2) Euler angle vector E.
%
*/
void	PRV2Euler232(double *q, double *e)
{
  double ep[5];
  
  PRV2EP(q,ep);
  EP2Euler232(ep,e);
  
  return;
}





/*
%
%	PRV2Euler312(Q,E) translates the principal rotation
%	vector Q into the (3-1-2) Euler angle vector E.
%
*/
void	PRV2Euler312(double *q, double *e)
{
  double ep[5];
  
  PRV2EP(q,ep);
  EP2Euler312(ep,e);
  
  return;
}





/*
%
%	PRV2Euler313(Q,E) translates the principal rotation
%	vector Q into the (3-1-3) Euler angle vector E.
%
*/
void	PRV2Euler313(double *q, double *e)
{
  double ep[5];
  
  PRV2EP(q,ep);
  EP2Euler313(ep,e);
  
  return;
}





/*
%
%	PRV2Euler321(Q,E) translates the principal rotation
%	vector Q into the (3-2-1) Euler angle vector E.
%
*/
void	PRV2Euler321(double *q, double *e)
{
  double ep[5];
  
  PRV2EP(q,ep);
  EP2Euler321(ep,e);
  
  return;
}






/*
%
%	PRV2Euler323(Q,E) translates the principal rotation
%	vector Q into the (3-2-3) Euler angle vector E.
%
*/
void	PRV2Euler323(double *q, double *e)
{
  double ep[5];
  
  PRV2EP(q,ep);
  EP2Euler323(ep,e);
  
  return;
}





/*
%
%	PRV2Gibbs(Q0,Q) translates the principal rotation vector Q1
%	into the Gibbs vector Q.
%
*/
void	PRV2Gibbs(double *q0, double *q)
{
  double q1[5], tp;
  
  PRV2elem(q0,q1);
  tp = tan(q1[1]/2.);
  q[1] = q1[2]*tp;
  q[2] = q1[3]*tp;
  q[3] = q1[4]*tp;
  
  return;
}





/*
%
%	PRV2MRP(Q0,Q) translates the principal rotation vector Q1
%	into the MRP vector Q.
%
*/
void	PRV2MRP(double *q0, double *q)
{
  double q1[5], tp;
  
  PRV2elem(q0,q1);
  tp = tan(q1[1]/4.);
  q[1] = q1[2]*tp;
  q[2] = q1[3]*tp;
  q[3] = q1[4]*tp;
  
  return;
}





/*
%
%	subEP(B1,B2,Q) provides the Euler parameter vector
%	which corresponds to relative rotation from B2
%	to B1.
%
*/
void	subEP(double *b1, double *b2, double *q)
{
  
  q[1] = b2[1]*b1[1]+b2[2]*b1[2]+b2[3]*b1[3]+b2[4]*b1[4];
  q[2] = -b2[2]*b1[1]+b2[1]*b1[2]+b2[4]*b1[3]-b2[3]*b1[4];
  q[3] = -b2[3]*b1[1]-b2[4]*b1[2]+b2[1]*b1[3]+b2[2]*b1[4];
  q[4] = -b2[4]*b1[1]+b2[3]*b1[2]-b2[2]*b1[3]+b2[1]*b1[4];
  
  return;
}




/*
%
%	subEuler121(E,E1,E2) computes the relative
%	(1-2-1) Euler angle vector from E1 to E.
%
*/
void	subEuler121(double *e, double *e1, double *e2)
{
  double cp, cp1, sp, sp1, cp2, dum;
  
  cp = cos(e[2]);
  cp1 = cos(e1[2]);
  sp = sin(e[2]);
  sp1 = sin(e1[2]);
  dum = e[1]-e1[1];

  e2[2] = acos(cp1*cp+sp1*sp*cos(dum));
  cp2 = cos(e2[2]);
  e2[1] = Picheck(-e1[3] + atan2(sp1*sp*sin(dum),cp2*cp1-cp));
  e2[3] = Picheck(e[3] - atan2(sp1*sp*sin(dum),cp1-cp*cp2));
  
  return;
}





/*
%
%	subEuler123(E,E1,E2) computes the relative
%	(1-2-3) Euler angle vector from E1 to E.
%
*/
void	subEuler123(double *e, double *e1, double *e2)
{
  double C[4][4], C1[4][4], C2[4][4];
  
  Euler1232C(e,C);
  Euler1232C(e1,C1);
  MdotMT(C,C1,C2);
  C2Euler123(C2,e2);
}




/*
%
%	subEuler131(E,E1,E2) computes the relative
%	(1-3-1) Euler angle vector from E1 to E.
%
*/
void	subEuler131(double *e, double *e1, double *e2)
{
  double cp, cp1, sp, sp1, dum, cp2;
  
  cp = cos(e[2]);
  cp1 = cos(e1[2]);
  sp = sin(e[2]);
  sp1 = sin(e1[2]);
  dum = e[1]-e1[1];

  e2[2] = acos(cp1*cp+sp1*sp*cos(dum));
  cp2 = cos(e2[2]);
  e2[1] = Picheck(-e1[3] + atan2(sp1*sp*sin(dum),cp2*cp1-cp));
  e2[3] = Picheck(e[3] - atan2(sp1*sp*sin(dum),cp1-cp*cp2));
  
  return;
}






/*
%
%	subEuler132(E,E1,E2) computes the relative
%	(1-3-2) Euler angle vector from E1 to E.
%
*/
void	subEuler132(double *e, double *e1, double *e2)
{
  double C[4][4], C1[4][4], C2[4][4];
  
  Euler1322C(e,C);
  Euler1322C(e1,C1);
  MdotMT(C,C1,C2);
  C2Euler132(C2,e2);
}





/*
%
%	subEuler212(E,E1,E2) computes the relative
%	(2-1-2) Euler angle vector from E1 to E.
%
*/
void	subEuler212(double *e, double *e1, double *e2)
{
  double cp, cp1, sp, sp1, dum, cp2;
  
  cp = cos(e[2]);
  cp1 = cos(e1[2]);
  sp = sin(e[2]);
  sp1 = sin(e1[2]);
  dum = e[1]-e1[1];

  e2[2] = acos(cp1*cp+sp1*sp*cos(dum));
  cp2 = cos(e2[2]);
  e2[1] = Picheck(-e1[3] + atan2(sp1*sp*sin(dum),cp2*cp1-cp));
  e2[3] = Picheck(e[3] - atan2(sp1*sp*sin(dum),cp1-cp*cp2));
  
  return;
}





/*
%
%	subEuler213(E,E1,E2) computes the relative
%	(2-1-3) Euler angle vector from E1 to E.
%
*/
void	subEuler213(double *e, double *e1, double *e2)
{
  double C[4][4], C1[4][4], C2[4][4];
  
  Euler2132C(e,C);
  Euler2132C(e1,C1);
  MdotMT(C,C1,C2);
  C2Euler213(C2,e2);
}








/*
%
%	subEuler231(E,E1,E2) computes the relative
%	(2-3-1) Euler angle vector from E1 to E.
%
*/
void	subEuler231(double *e, double *e1, double *e2)
{
  double C[4][4], C1[4][4], C2[4][4];
  
  Euler2312C(e,C);
  Euler2312C(e1,C1);
  MdotMT(C,C1,C2);
  C2Euler231(C2,e2);
}




/*
%
%	subEuler232(E,E1,E2) computes the relative
%	(2-3-2) Euler angle vector from E1 to E.
%
*/
void	subEuler232(double *e, double *e1, double *e2)
{
  double cp, cp1, sp, sp1, dum, cp2;
  
  cp = cos(e[2]);
  cp1 = cos(e1[2]);
  sp = sin(e[2]);
  sp1 = sin(e1[2]);
  dum = e[1]-e1[1];

  e2[2] = acos(cp1*cp+sp1*sp*cos(dum));
  cp2 = cos(e2[2]);
  e2[1] = Picheck(-e1[3] + atan2(sp1*sp*sin(dum),cp2*cp1-cp));
  e2[3] = Picheck(e[3] - atan2(sp1*sp*sin(dum),cp1-cp*cp2));
  
  return;
}





/*
%
%	subEuler312(E,E1,E2) computes the relative
%	(3-1-2) Euler angle vector from E1 to E.
%
*/
void	subEuler312(double *e, double *e1, double *e2)
{
  double C[4][4], C1[4][4], C2[4][4];
  
  Euler3122C(e,C);
  Euler3122C(e1,C1);
  MdotMT(C,C1,C2);
  C2Euler312(C2,e2);
}








/*
%
%	subEuler313(E,E1,E2) computes the relative
%	(3-1-3) Euler angle vector from E1 to E.
%
*/
void	subEuler313(double *e, double *e1, double *e2)
{
  double cp, cp1, sp, sp1, dum, cp2;
  
  cp = cos(e[2]);
  cp1 = cos(e1[2]);
  sp = sin(e[2]);
  sp1 = sin(e1[2]);
  dum = e[1]-e1[1];

  e2[2] = acos(cp1*cp+sp1*sp*cos(dum));
  cp2 = cos(e2[2]);
  e2[1] = Picheck(-e1[3] + atan2(sp1*sp*sin(dum),cp2*cp1-cp));
  e2[3] = Picheck(e[3] - atan2(sp1*sp*sin(dum),cp1-cp*cp2));
  
  return;
}






/*
%
%	subEuler321(E,E1,E2) computes the relative
%	(3-2-1) Euler angle vector from E1 to E.
%
*/
void	subEuler321(double *e, double *e1, double *e2)
{
  double C[4][4], C1[4][4], C2[4][4];
  
  Euler3212C(e,C);
  Euler3212C(e1,C1);
  MdotMT(C,C1,C2);
  C2Euler321(C2,e2);
}





/*
%
%	subEuler323(E,E1,E2) computes the relative
%	(3-2-3) Euler angle vector from E1 to E.
%
*/
void	subEuler323(double *e, double *e1, double *e2)
{
  double cp, cp1, sp, sp1, dum, cp2;
  
  cp = cos(e[2]);
  cp1 = cos(e1[2]);
  sp = sin(e[2]);
  sp1 = sin(e1[2]);
  dum = e[1]-e1[1];

  e2[2] = acos(cp1*cp+sp1*sp*cos(dum));
  cp2 = cos(e2[2]);
  e2[1] = Picheck(-e1[3] + atan2(sp1*sp*sin(dum),cp2*cp1-cp));
  e2[3] = Picheck(e[3] - atan2(sp1*sp*sin(dum),cp1-cp*cp2));
  
  return;
}






/*
%
%	subGibbs(Q1,Q2,Q) provides the Gibbs vector
%	which corresponds to relative rotation from Q2
%	to Q1.
%
*/
void	subGibbs(double *q1, double *q2, double *q)
{
  double d1[4];
  
  cross(q1,q2,d1);
  add(q1,d1,q);
  sub(q,q2,q);
  mult(1./(1.+dot(q1,q2)),q,q);
  
  return;
}





/*
%
%	subMRP(Q1,Q2,Q) provides the MRP vector
%	which corresponds to relative rotation from Q2
%	to Q1.
%
*/
void	subMRP(double *q1, double *q2, double *q)
{
  double d1[4];
  
  cross(q1,q2,d1);
  mult(2.,d1,q);
  mult(1.-dot(q2,q2),q1,d1);
  add(q,d1,q);
  mult(1.-dot(q1,q1),q2,d1);
  sub(q,d1,q);
  mult(1./(1.+dot(q1,q1)*dot(q2,q2)+2.*dot(q1,q2)),q,q);
  
  return;
}






/*
%
%	subPRV(Q1,Q2,Q) provides the prinipal rotation vector
%	which corresponds to relative principal rotation from Q2
%	to Q1.
%
*/
void	subPRV(double *q10, double *q20, double *q)
{
  double q1[5], q2[5], cp1, cp2, sp1, sp2, e1[4], e2[4], p, sp;
  
  PRV2elem(q10,q1);
  PRV2elem(q20,q2);
  cp1 = cos(q1[1]/2.);
  cp2 = cos(q2[1]/2.);
  sp1 = sin(q1[1]/2.);
  sp2 = sin(q2[1]/2.);
  equal(&(q1[1]), e1);
  equal(&(q2[1]), e2);
  
  p = 2.*acos(cp1*cp2+sp1*sp2*dot(e1,e2));
  sp = sin(p/2.);
  
  cross(e1,e2,q1);
  mult(sp1*sp2,q1,q);
  mult(cp2*sp1,e1,q1);
  add(q1,q,q);
  mult(cp1*sp2,e2,q1);
  sub(q,q1,q);
  mult(p/sp,q,q);
    
  return;
}


/*
%
%	Mi(theta, a, C) returns the rotation matrix corresponding
%	to a single axis rotation about axis a by the angle theta
%
*/
void  Mi(double theta, int a, double C[4][4])
{
  double c,s;
  
  c = cos(theta);
  s = sin(theta);
  
  switch (a) {
    case 1:
      C[1][1] = 1.;   C[1][2] = 0.;   C[1][3] = 0.;
      C[2][1] = 0.;   C[2][2] =  c;   C[2][3] =  s; 
      C[3][1] = 0.;   C[3][2] = -s;   C[3][3] =  c;
      break;
      
    case 2:
      C[1][1] =  c;   C[1][2] = 0.;   C[1][3] = -s;
      C[2][1] = 0.;   C[2][2] = 1.;   C[2][3] = 0.; 
      C[3][1] =  s;   C[3][2] = 0.;   C[3][3] =  c;
      break;
    
    case 3:
      C[1][1] =  c;   C[1][2] =  s;   C[1][3] = 0.;
      C[2][1] = -s;   C[2][2] =  c;   C[2][3] = 0.; 
      C[3][1] = 0.;   C[3][2] = 0.;   C[3][3] = 1.;
      break;
      
    default:
      printf("Mi() error: incorrect axis %d selected.\n", a);
  }
  
  return;
}

