(* :Title: RigidBodyKinematics *)

(* :Author: Hanspeter Schaub *)

(* :Summary:
This package provides various rigid body kinematic routines.  Functions
are included that translate between various attitude parameters, computes 
the successive rotation combination and the differential, kinematic 
equation.  

	Abbreviations:
		C		- Direction Cosine Matrix
		PRV		- Prinicpal Rotation Vector
		MRP		- Modified Rodrigues Parameter vector
		Gibbs		- Classical Rodrigues Parameter vector or Gibbs vector
		Eulerijk	- (i-j-k) Euler angles
*)


(* :Package Version: 1.0.1 *)

(* :Copyright: Freeware *)

(* :History: 
	Version 1.0 by Hanspeter Schaub (Texas A&M University), February 1999.
	Version 1.0.1 by Hanspeter Schaub (Virginia Tech), June 2005.
*)

(* :Keywords: rigid body kinematics, differential equation *)

(* :Mathematica Version: 3.0 *)

BeginPackage["RigidBodyKinematics`"]

(* Usage messages *)

addEP::usage =
	"Q = addEP[B1,B2] provides the Euler parameter vector
	which corresponds to performing to successive
	rotations B1 and B2."

addEuler121::usage = 
	"Q = addEuler121[E1,E2] computes the overall (1-2-1) Euler
	angle vector corresponding to two successive
	(1-2-1) rotations E1 and E2."

addEuler123::usage = 
	"Q = addEuler123[E1,E2] computes the overall (1-2-3) Euler
	angle vector corresponding to two successive
	(1-2-3) rotations E1 and E2."

addEuler131::usage = 
	"Q = addEuler131[E1,E2] computes the overall (1-3-1) Euler
	angle vector corresponding to two successive
	(1-3-1) rotations E1 and E2."

addEuler132::usage = 
	"Q = addEuler132[E1,E2] computes the overall (1-3-2) Euler
	angle vector corresponding to two successive
	(1-3-2) rotations E1 and E2."
	
addEuler212::usage = 
	"Q = addEuler212[E1,E2] computes the overall (2-1-2) Euler
	angle vector corresponding to two successive
	(2-1-2) rotations E1 and E2."

addEuler213::usage = 
	"Q = addEuler213[E1,E2] computes the overall (2-1-3) Euler
	angle vector corresponding to two successive
	(2-1-3) rotations E1 and E2."
	
addEuler231::usage = 
	"Q = addEuler231[E1,E2] computes the overall (2-3-1) Euler
	angle vector corresponding to two successive
	(2-3-1) rotations E1 and E2."

addEuler232::usage = 
	"Q = addEuler232[E1,E2] computes the overall (2-3-2) Euler
	angle vector corresponding to two successive
	(2-3-2) rotations E1 and E2."

addEuler312::usage = 
	"Q = addEuler312[E1,E2] computes the overall (3-1-2) Euler
	angle vector corresponding to two successive
	(3-1-2) rotations E1 and E2."

addEuler313::usage = 
	"Q = addEuler313[E1,E2] computes the overall (3-1-3) Euler
	angle vector corresponding to two successive
	(3-1-3) rotations E1 and E2."

addEuler321::usage = 
	"Q = addEuler321[E1,E2] computes the overall (3-2-1) Euler
	angle vector corresponding to two successive
	(3-2-1) rotations E1 and E2."

addEuler323::usage = 
	"Q = addEuler323[E1,E2] computes the overall (3-2-3) Euler
	angle vector corresponding to two successive
	(3-2-3) rotations E1 and E2."

addGibbs::usage = 
	"Q = addGibbs[Q1,Q2] provides the Gibbs vector
	which corresponds to performing to successive
	rotations Q1 and Q2."

addMRP::usage = 
	"Q = addMRP[Q1,Q2] provides the MRP vector
	which corresponds to performing to successive
	rotations Q1 and Q2."

addPRV::usage = 
	"Q = addPRV[Q1,Q2] provides the principal rotation vector
	which corresponds to performing to successive
	prinicipal rotations Q1 and Q2."

BinvEP::usage = 
	"B = BinvEP[Q] returns the 3x4 matrix which relates 
	the derivative of Euler parameter vector Q to the 
	body angular velocity vector w.  
	
		w = 2 [B(Q)]^(-1) dQ/dt"

BinvEuler121::usage = 
	"B = BinvEuler121[Q] returns the 3x3 matrix which relates 
	the derivative of the (1-2-1) Euler angle vector Q to the 
	body angular velocity vector w.  	
	
		w = [B(Q)]^(-1) dQ/dt"

BinvEuler123::usage = 
	"B = BinvEuler123[Q] returns the 3x3 matrix which relates 
	the derivative of the (1-2-3) Euler angle vector Q to the 
	body angular velocity vector w.  
	
		w = [B(Q)]^(-1) dQ/dt"

BinvEuler131::usage =
	"B = BinvEuler131[Q] returns the 3x3 matrix which relates 
	the derivative of the (1-3-1) Euler angle vector Q to the 
	body angular velocity vector w.  
	
		w = [B(Q)]^(-1) dQ/dt"

BinvEuler132::usage = 
	"B = BinvEuler132[Q] returns the 3x3 matrix which relates 
	the derivative of the (1-3-2) Euler angle vector Q to the 
	body angular velocity vector w.  
	
		w = [B(Q)]^(-1) dQ/dt"

BinvEuler212::usage = 
	"B = BinvEuler212[Q] returns the 3x3 matrix which relates 
	the derivative of the (2-1-2) Euler angle vector Q to the 
	body angular velocity vector w.  
	
		w = [B(Q)]^(-1) dQ/dt"

BinvEuler213::usage = 
	"B = BinvEuler213[Q] returns the 3x3 matrix which relates 
	the derivative of the (2-1-3) Euler angle vector Q to the 
	body angular velocity vector w.  
	
		w = [B(Q)]^(-1) dQ/dt"

BinvEuler231::usage = 
	"B = BinvEuler231[Q] returns the 3x3 matrix which relates 
	the derivative of the (2-3-1) Euler angle vector Q to the 
	body angular velocity vector w.  
	
		w = [B(Q)]^(-1) dQ/dt"

BinvEuler232::usage =
	"B = BinvEuler232[Q] returns the 3x3 matrix which relates 
	the derivative of the (2-3-2) Euler angle vector Q to the 
	body angular velocity vector w.  
	
		w = [B(Q)]^(-1) dQ/dt"

BinvEuler312::usage = 
	"B = BinvEuler312[Q] returns the 3x3 matrix which relates 
	the derivative of the (3-1-2) Euler angle vector Q to the 
	body angular velocity vector w.  
	
		w = [B(Q)]^(-1) dQ/dt"

BinvEuler313::usage = 
	"B = BinvEuler313[Q] returns the 3x3 matrix which relates 
	the derivative of the (3-1-3) Euler angle vector Q to the 
	body angular velocity vector w.  
	
		w = [B(Q)]^(-1) dQ/dt"

BinvEuler321::usage = 
	"B = BinvEuler321[Q] returns the 3x3 matrix which relates 
	the derivative of the (3-2-1) Euler angle vector Q to the 
	body angular velocity vector w.  
	
		w = [B(Q)]^(-1) dQ/dt"

BinvEuler323::usage = 
	"B = BinvEuler323[Q] returns the 3x3 matrix which relates 
	the derivative of the (3-2-3) Euler angle vector Q to the 
	body angular velocity vector w.  
	
		w = [B(Q)]^(-1) dQ/dt"

BinvGibbs::usage = 
	"B = BinvGibbs[Q] returns the 3x3 matrix which relates 
	the derivative of Gibbs vector Q to the 
	body angular velocity vector w.  
	
		w = 2 [B(Q)]^(-1) dQ/dt"

BinvMRP::usage = 
	"B = BinvMRP[Q] returns the 3x3 matrix which relates 
	the derivative of MRP vector Q to the 
	body angular velocity vector w.  
	
		w = 4 [B(Q)]^(-1) dQ/dt"

BinvPRV::usage = 
	"B = BinvPRV[Q] returns the 3x3 matrix which relates 
	the derivative of principal rotation vector Q to the 
	body angular velocity vector w.  
	
		w = [B(Q)]^(-1) dQ/dt"

BmatEP::usage = 
	"B = BmatEP[Q] returns the 4x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	Euler parameter vector Q.  
	
		dQ/dt = 1/2 [B(Q)] w"

BmatEuler121::usage = 
	"B = BmatEuler121[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	(1-2-1) Euler angle vector Q.  
	
		dQ/dt = [B(Q)] w"

BmatEuler123::usage = 
	"B = BmatEuler123[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	(1-2-3) Euler angle vector Q.  
	
		dQ/dt = [B(Q)] w"

BmatEuler131::usage = 
	"B = BmatEuler131[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	(1-3-1) Euler angle vector Q.  
	
		dQ/dt = [B(Q)] w"

BmatEuler132::usage = 
	"B = BmatEuler132[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	(1-3-2) Euler angle vector Q.  
	
		dQ/dt = [B(Q)] w"

BmatEuler212::usage = 
	"B = BmatEuler212[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	(2-1-2) Euler angle vector Q.  
	
		dQ/dt = [B(Q)] w"

BmatEuler213::usage = 
	"B = BmatEuler213[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	(2-1-3) Euler angle vector Q.  
	
		dQ/dt = [B(Q)] w"

BmatEuler231::usage = 
	"B = BmatEuler231[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	(2-3-1) Euler angle vector Q.  
	
		dQ/dt = [B(Q)] w"

BmatEuler232::usage = 
	"B = BmatEuler232[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	(2-3-2) Euler angle vector Q.  
	
		dQ/dt = [B(Q)] w"

BmatEuler312::usage = 
	"B = BmatEuler312[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	(3-1-2) Euler angle vector Q.  
	
		dQ/dt = [B(Q)] w"

BmatEuler313::usage = 
	"B = BmatEuler313[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	(3-1-3) Euler angle vector Q.  
	
		dQ/dt = [B(Q)] w"

BmatEuler321::usage = 
	"B = BmatEuler321[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	(3-2-1) Euler angle vector Q.  
	
		dQ/dt = [B(Q)] w"

BmatEuler323::usage = 
	"B = BmatEuler323[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	(3-2-3) Euler angle vector Q.  
	
		dQ/dt = [B(Q)] w"

BmatGibbs::usage = 	
	"B = BmatGibbs[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	Gibbs vector Q.  
	
		dQ/dt = 1/2 [B(Q)] w"

BmatMRP::usage = 
	"B = BmatMRP[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	MRP vector Q.  
	
		dQ/dt = 1/4 [B(Q)] w"

BmatPRV::usage = 
	"B = BmatPRV[Q] returns the 3x3 matrix which relates the 
	body angular velocity vector w to the derivative of
	principal rotation vector Q.  
	
		dQ/dt = [B(Q)] w"

C2EP::usage = 
	"Q = C2EP[C] translates the 3x3 direction cosine matrix
	C into the corresponding 4x1 Euler parameter vector Q,
	where the first component of Q is the non-dimensional
	Euler parameter Beta_0 >= 0.  Transformation is done
	using the Stanley method."

C2Euler121::usage = 
	"Q = C2Euler121[C] translates the 3x3 direction cosine matrix
	C into the corresponding (1-2-1) Euler angle set."

C2Euler123::usage = 
	"Q = C2Euler123[C] translates the 3x3 direction cosine matrix
	C into the corresponding (1-2-3) Euler angle set."

C2Euler131::usage = 
	"Q = C2Euler131[C] translates the 3x3 direction cosine matrix
	C into the corresponding (1-3-1) Euler angle set."

C2Euler132::usage =
	"Q = C2Euler132[C] translates the 3x3 direction cosine matrix
	C into the corresponding (1-3-2) Euler angle set."

C2Euler212::usage = 
	"Q = C2Euler212[C] translates the 3x3 direction cosine matrix
	C into the corresponding (2-1-2) Euler angle set."

C2Euler213::usage = 
	"Q = C2Euler213[C]) translates the 3x3 direction cosine matrix
	C into the corresponding (2-1-3) Euler angle set."

C2Euler231::usage = 
	"Q = C2Euler231[C] translates the 3x3 direction cosine matrix
	C into the corresponding (2-3-1) Euler angle set."

C2Euler232::usage = 
	"Q = C2Euler232[C] translates the 3x3 direction cosine matrix
	C into the corresponding (2-3-2) Euler angle set."

C2Euler312::usage = 
	"Q = C2Euler312[C] translates the 3x3 direction cosine matrix
	C into the corresponding (3-1-2) Euler angle set."

C2Euler313::usage = 
	"Q = C2Euler313[C] translates the 3x3 direction cosine matrix
	C into the corresponding (3-1-3) Euler angle set."

C2Euler321::usage = 
	"Q = C2Euler321[C] translates the 3x3 direction cosine matrix
	C into the corresponding (3-2-1) Euler angle set."

C2Euler323::usage = 
	"Q = C2Euler323[C] translates the 3x3 direction cosine matrix
	C into the corresponding (3-2-3) Euler angle set."

C2Gibbs::usage = 
	"Q = C2Gibbs[C] translates the 3x3 direction cosine matrix
	C into the corresponding 3x1 Gibbs vector Q."

C2MRP::usage = 
	"Q = C2MRP[C] translates the 3x3 direction cosine matrix
	C into the corresponding 3x1 MRP vector Q where the 
	MRP vector is chosen such that |Q| <= 1."

C2PRV::usage = 	
	"Q = C2PRV[C] translates the 3x3 direction cosine matrix
	C into the corresponding 3x1 principal rotation vector Q,
	where the first component of Q is the principal rotation angle
	phi (0<= phi <= Pi)"

dEP::usage = 	
	"dq = dEP(Q,W) returns the Euler parameter derivative
	for a given Euler parameter vector Q and body
	angular velocity vector w.

	dQ/dt = 1/2 [B(Q)] w"
	
dEuler121::usage = 	
	"dq = dEuler121[Q,W] returns the (1-2-1) Euler angle derivative
	vector for a given (1-2-1) Euler angle vector Q and body
	angular velocity vector w.

	dQ/dt =  [B(Q)] w"
	
dEuler123::usage = 	
	"dq = dEuler123[Q,W] returns the (1-2-3) Euler angle derivative
	vector for a given (1-2-3) Euler angle vector Q and body
	angular velocity vector w.

	dQ/dt =  [B(Q)] w"

dEuler131::usage = 	
	"dq = dEuler131[Q,W] returns the (1-3-1) Euler angle derivative
	vector for a given (1-3-1) Euler angle vector Q and body
	angular velocity vector w.

	dQ/dt =  [B(Q)] w"

dEuler132::usage = 	
	"dq = dEuler132[Q,W] returns the (1-3-2) Euler angle derivative
	vector for a given (1-3-2) Euler angle vector Q and body
	angular velocity vector w.

	dQ/dt =  [B(Q)] w"

dEuler232::usage = 	
	"dq = dEuler232[Q,W] returns the (2-3-2) Euler angle derivative
	vector for a given (2-3-2) Euler angle vector Q and body
	angular velocity vector w.

	dQ/dt =  [B(Q)] w"

dEuler231::usage = 	
	"dq = dEuler231[Q,W] returns the (2-3-1) Euler angle derivative
	vector for a given (2-3-1) Euler angle vector Q and body
	angular velocity vector w.

	dQ/dt =  [B(Q)] w"

dEuler212::usage = 	
	"dq = dEuler212[Q,W] returns the (2-1-2) Euler angle derivative
	vector for a given (2-1-2) Euler angle vector Q and body
	angular velocity vector w.

	dQ/dt =  [B(Q)] w"

dEuler213::usage = 	
	"dq = dEuler213[Q,W] returns the (2-1-3) Euler angle derivative
	vector for a given (2-1-3) Euler angle vector Q and body
	angular velocity vector w.

	dQ/dt =  [B(Q)] w"

dEuler313::usage = 	
	"dq = dEuler313[Q,W] returns the (3-1-3) Euler angle derivative
	vector for a given (3-1-3) Euler angle vector Q and body
	angular velocity vector w.

	dQ/dt =  [B(Q)] w"

dEuler312::usage = 	
	"dq = dEuler312[Q,W] returns the (3-1-2) Euler angle derivative
	vector for a given (3-1-2) Euler angle vector Q and body
	angular velocity vector w.

	dQ/dt =  [B(Q)] w"

dEuler323::usage = 	
	"dq = dEuler323[Q,W] returns the (3-2-3) Euler angle derivative
	vector for a given (3-2-3) Euler angle vector Q and body
	angular velocity vector w.

	dQ/dt =  [B(Q)] w"

dEuler321::usage = 	
	"dq = dEuler321[Q,W] returns the (3-2-1) Euler angle derivative
	vector for a given (3-2-1) Euler angle vector Q and body
	angular velocity vector w.

	dQ/dt =  [B(Q)] w"

dGibbs::usage = 
	"dq = dGibbs[Q,W] returns the Gibbs derivative
	for a given Gibbs vector Q and body
	angular velocity vector w.

	dQ/dt = 1/2 [B(Q)] w"

dMRP::usage = 
	"dq = dMRP[Q,W] returns the MRP derivative
	for a given MRP vector Q and body
	angular velocity vector w.

	dQ/dt = 1/4 [B(Q)] w"

dPRV::usage = 
	"dq = dPRV[Q,W] returns the PRV derivative
	for a given PRV vector Q and body
	angular velocity vector w.

	dQ/dt =  [B(Q)] w"

elem2PRV::usage = 
	"Q = elem2PRV[R] translates a prinicpal rotation 
	element set R into the corresponding principal 
	rotation vector Q."

EP2C::usage = 
	"C = EP2C[Q] returns the direction cosine 
	matrix in terms of the 4x1 Euler parameter vector
	Q.  The first element is the non-dimensional Euler
	parameter, while the remain three elements form 
	the Eulerparameter vector."

EP2Euler121::usage = 
	"E = EP2Euler121[Q] translates the Euler parameter
	vector Q into the corresponding (1-2-1) Euler angle
	vector E."

EP2Euler123::usage = 
	"E = EP2Euler123[Q] translates the Euler parameter vector
	Q into the corresponding (1-2-3) Euler angle set."

EP2Euler131::usage = 
	"E = EP2Euler131[Q] translates the Euler parameter
	vector Q into the corresponding (1-3-1) Euler angle
	vector E."

EP2Euler132::usage = 
	"E = EP2Euler132[Q] translates the Euler parameter vector
	Q into the corresponding (1-3-2) Euler angle set."

EP2Euler212::usage = 
	"E = EP2Euler212[Q] translates the Euler parameter
	vector Q into the corresponding (2-1-2) Euler angle
	vector E."

EP2Euler213::usage = 
	"E = EP2Euler213[Q] translates the Euler parameter vector
	Q into the corresponding (2-1-3) Euler angle set."

EP2Euler231::usage = 
	"E = EP2Euler231[Q] translates the Euler parameter vector
	Q into the corresponding (2-3-1) Euler angle set."

EP2Euler232::usage = 	
	"E = EP2Euler232[Q] translates the Euler parameter
	vector Q into the corresponding (2-3-2) Euler angle
	vector E."

EP2Euler312::usage = 
	"E = EP2Euler312[Q] translates the Euler parameter vector
	Q into the corresponding (3-1-2) Euler angle set."

EP2Euler313::usage = 
	"E = EP2Euler313[Q] translates the Euler parameter
	vector Q into the corresponding (3-1-3) Euler angle
	vector E."

EP2Euler321::usage = 
	"E = EP2Euler321[Q] translates the Euler parameter vector
	Q into the corresponding (3-2-1) Euler angle set."

EP2Euler323::usage = 
	"E = EP2Euler323[Q] translates the Euler parameter
	vector Q into the corresponding (3-2-3) Euler angle
	vector E."

EP2Gibbs::usage = 
	"Q = EP2Gibbs[Q1] translates the Euler parameter vector Q1
	into the Gibbs vector Q."

EP2MRP::usage = 
	"Q = EP2MRP[Q1] translates the Euler parameter vector Q1
	into the MRP vector Q."

EP2PRV::usage = 
	"Q = EP2PRV[Q1] translates the Euler parameter vector Q1
	into the principal rotation vector Q."

Euler1::usage = 
	"C = Euler[x] Returns the elementary rotation matrix about the
	first body axis."

Euler2::usage = 
	"C = Euler2[x] Returns the elementary rotation matrix about the
	second body axis."

Euler3::usage = 
	"C = Euler3(x) Returns the elementary rotation matrix about the
	third body axis."

Euler1212C::usage = 
	"C = Euler1212C[Q] returns the direction cosine 
	matrix in terms of the 1-2-1 Euler angles.  
	Input Q must be a 3x1 vector of Euler angles."

Euler1232C::usage = 
	"C = Euler1232C[Q] returns the direction cosine 
	matrix in terms of the 1-2-3 Euler angles.  
	Input Q must be a 3x1 vector of Euler angles."

Euler1312C::usage = 
	"C = Euler1312C[Q] returns the direction cosine 
	matrix in terms of the 1-3-1 Euler angles.  
	Input Q must be a 3x1 vector of Euler angles."

Euler1322C::usage =
	"C = Euler1322C[Q] returns the direction cosine 
	matrix in terms of the 1-3-2 Euler angles.  
	Input Q must be a 3x1 vector of Euler angles."

Euler2122C::usage =
	"C = Euler2122C[Q] returns the direction cosine 
	matrix in terms of the 2-1-2 Euler angles.  
	Input Q must be a 3x1 vector of Euler angles."

Euler2132C::usage = 
	"C = Euler2132C[Q] returns the direction cosine 
	matrix in terms of the 2-1-3 Euler angles.  
	Input Q must be a 3x1 vector of Euler angles."

Euler2312C::usage = 
	"C = Euler2312C[Q] returns the direction cosine 
	matrix in terms of the 2-3-1 Euler angles.  
	Input Q must be a 3x1 vector of Euler angles."

Euler2322C::usage = 
	"C = Euler2322C[Q] returns the direction cosine 
	matrix in terms of the 2-3-2 Euler angles.  
	Input Q must be a 3x1 vector of Euler angles."

Euler3122C::usage = 
	"C = Euler3122C[Q] returns the direction cosine 
	matrix in terms of the 3-1-2 Euler angles.  
	Input Q must be a 3x1 vector of Euler angles."

Euler3132C::usage = 
	"C = Euler3132C[Q] returns the direction cosine 
	matrix in terms of the 3-1-3 Euler angles.  
	Input Q must be a 3x1 vector of Euler angles."

Euler3212C::usage = 
	"C = Euler3212C[Q] returns the direction cosine 
	matrix in terms of the 3-2-1 Euler angles.  
	Input Q must be a 3x1 vector of Euler angles."

Euler3232C::usage =
	"C = Euler3232C[Q] returns the direction cosine 
	matrix in terms of the 3-2-3 Euler angles.  
	Input Q must be a 3x1 vector of Euler angles."

Euler1212EP::usage = 
	"Q = Euler1212EP[E] translates the 121 Euler angle
	vector E into the Euler parameter vector Q."

Euler1232EP::usage = 
	"Q = Euler1232EP[E] translates the 123 Euler angle
	vector E into the Euler parameter vector Q."

Euler1312EP::usage = 
	"Q = Euler1312EP[E] translates the 131 Euler angle
	vector E into the Euler parameter vector Q."

Euler1322EP::usage = 
	"Q = Euler1322EP[E] translates the 132 Euler angle
	vector E into the Euler parameter vector Q."

Euler2122EP::usage = 
	"Q = Euler2122EP[E] translates the 212 Euler angle
	vector E into the Euler parameter vector Q."

Euler2132EP::usage = 
	"Q = Euler2132EP[E] translates the 213 Euler angle
	vector E into the Euler parameter vector Q."

Euler2312EP::usage = 
	"Q = Euler2312EP[E] translates the 231 Euler angle
	vector E into the Euler parameter vector Q."

Euler2322EP::usage = 
	"Q = Euler2322EP[E] translates the 232 Euler angle
	vector E into the Euler parameter vector Q."

Euler3122EP::usage = 
	"Q = Euler3122EP[E] translates the 312 Euler angle
	vector E into the Euler parameter vector Q."

Euler3132EP::usage =
	"Q = Euler3132EP[E] translates the 313 Euler angle
	vector E into the Euler parameter vector Q."

Euler3212EP::usage =
	"Q = Euler2132EP[E] translates the 321 Euler angle
	vector E into the Euler parameter vector Q."

Euler3232EP::usage =
	"Q = Euler3232EP(E) translates the 323 Euler angle
	vector E into the Euler parameter vector Q."

Euler1212Gibbs::usage =
	"Q = Euler1212Gibbs[E] translates the (1-2-1) Euler
	angle vector E into the Gibbs vector Q."

Euler1232Gibbs::usage =
	"Q = Euler1232Gibbs[E] translates the (1-2-3) Euler
	angle vector E into the Gibbs vector Q."

Euler1312Gibbs::usage =
	"Q = Euler1312Gibbs[E] translates the (1-3-1) Euler
	angle vector E into the Gibbs vector Q."

Euler1322Gibbs::usage =
	"Q = Euler1322Gibbs[E] translates the (1-3-2) Euler
	angle vector E into the Gibbs vector Q."

Euler2122Gibbs::usage =
	"Q = Euler2122Gibbs[E] translates the (2-1-2) Euler
	angle vector E into the Gibbs vector Q."

Euler2132Gibbs::usage =
	"Q = Euler2132Gibbs[E] translates the (2-1-3) Euler
	angle vector E into the Gibbs vector Q."

Euler2312Gibbs::usage =
	"Q = Euler2312Gibbs[E] translates the (2-3-1) Euler
	angle vector E into the Gibbs vector Q."

Euler2322Gibbs::usage =
	"Q = Euler2322Gibbs[E] translates the (2-3-2) Euler
	angle vector E into the Gibbs vector Q."

Euler3132Gibbs::usage =
	"Q = Euler3132Gibbs[E] translates the (3-1-3) Euler
	angle vector E into the Gibbs vector Q."

Euler3122Gibbs::usage =
	"Q = Euler3122Gibbs[E] translates the (3-1-2) Euler
	angle vector E into the Gibbs vector Q."

Euler3232Gibbs::usage =
	"Q = Euler3232Gibbs[E] translates the (3-2-3) Euler
	angle vector E into the Gibbs vector Q."

Euler3212Gibbs::usage =
	"Q = Euler3212Gibbs[E] translates the (3-2-1) Euler
	angle vector E into the Gibbs vector Q."

Euler1212MRP::usage = 
	"Q = Euler1212MRP[E] translates the (1-2-1) Euler
	angle vector E into the MRP vector Q."

Euler1232MRP::usage = 
	"Q = Euler1232MRP[E] translates the (1-2-3) Euler
	angle vector E into the MRP vector Q."

Euler1312MRP::usage = 
	"Q = Euler1312MRP[E] translates the (1-3-1) Euler
	angle vector E into the MRP vector Q."

Euler1322MRP::usage = 
	"Q = Euler1322MRP[E] translates the (1-3-2) Euler
	angle vector E into the MRP vector Q."

Euler2122MRP::usage = 
	"Q = Euler2122MRP[E] translates the (2-1-2) Euler
	angle vector E into the MRP vector Q."

Euler2132MRP::usage = 
	"Q = Euler2132MRP[E] translates the (2-1-3) Euler
	angle vector E into the MRP vector Q."

Euler2322MRP::usage = 
	"Q = Euler2322MRP[E] translates the (2-3-2) Euler
	angle vector E into the MRP vector Q."

Euler2312MRP::usage = 
	"Q = Euler2312MRP[E] translates the (2-3-1) Euler
	angle vector E into the MRP vector Q."

Euler3132MRP::usage = 
	"Q = Euler3132MRP[E] translates the (3-1-3) Euler
	angle vector E into the MRP vector Q."

Euler3122MRP::usage = 
	"Q = Euler3122MRP[E] translates the (3-1-2) Euler
	angle vector E into the MRP vector Q."

Euler3232MRP::usage = 
	"Q = Euler3232MRP[E] translates the (3-2-3) Euler
	angle vector E into the MRP vector Q."

Euler3212MRP::usage = 
	"Q = Euler3212MRP[E] translates the (3-2-1) Euler
	angle vector E into the MRP vector Q."

Euler1212PRV::usage =
	"Q = Euler1212PRV[E] translates the (1-2-1) Euler
	angle vector E into the PRV vector Q."

Euler1232PRV::usage =
	"Q = Euler1232PRV[E] translates the (1-2-3) Euler
	angle vector E into the PRV vector Q."

Euler1312PRV::usage =
	"Q = Euler1312PRV[E] translates the (1-3-1) Euler
	angle vector E into the PRV vector Q."

Euler1322PRV::usage =
	"Q = Euler1322PRV[E] translates the (1-3-2) Euler
	angle vector E into the PRV vector Q."

Euler2122PRV::usage =
	"Q = Euler2122PRV[E] translates the (2-1-2) Euler
	angle vector E into the PRV vector Q."

Euler2132PRV::usage =
	"Q = Euler2132PRV[E] translates the (2-1-3) Euler
	angle vector E into the PRV vector Q."

Euler2312PRV::usage =
	"Q = Euler2312PRV[E] translates the (2-3-1) Euler
	angle vector E into the PRV vector Q."

Euler2322PRV::usage =
	"Q = Euler2322PRV[E] translates the (2-3-2) Euler
	angle vector E into the PRV vector Q."

Euler3132PRV::usage =
	"Q = Euler3132PRV[E] translates the (3-1-3) Euler
	angle vector E into the PRV vector Q."

Euler3122PRV::usage =
	"Q = Euler3122PRV[E] translates the (3-1-2) Euler
	angle vector E into the PRV vector Q."

Euler3232PRV::usage =
	"Q = Euler3232PRV[E] translates the (3-2-3) Euler
	angle vector E into the PRV vector Q."

Euler3212PRV::usage =
	"Q = Euler3212PRV[E] translates the (3-2-1) Euler
	angle vector E into the PRV vector Q."

Gibbs2C::usage = 
	"C = Gibbs2C[Q] returns the direction cosine 
	matrix in terms of the 3x1 Gibbs vector Q."

Gibbs2EP::usage = 
	"Q = Gibbs2EP[Q1] translates the Gibbs vector Q1
	into the Euler parameter vector Q."

Gibbs2Euler121::usage = 
	"E = Gibbs2Euler121[Q] translates the Gibbs
	 vector Q into the (1-2-1) Euler angle vector E."
	
Gibbs2Euler123::usage = 
	"E = Gibbs2Euler123[Q] translates the Gibbs
	 vector Q into the (1-2-3) Euler angle vector E."
	
Gibbs2Euler131::usage = 
	"E = Gibbs2Euler131[Q] translates the Gibbs
	 vector Q into the (1-3-1) Euler angle vector E."
	
Gibbs2Euler132::usage = 
	"E = Gibbs2Euler132[Q] translates the Gibbs
	 vector Q into the (1-3-2) Euler angle vector E."
	
Gibbs2Euler321::usage = 
	"E = Gibbs2Euler321[Q] translates the Gibbs
	 vector Q into the (3-2-1) Euler angle vector E."
	
Gibbs2Euler323::usage = 
	"E = Gibbs2Euler323[Q] translates the Gibbs
	 vector Q into the (3-2-3) Euler angle vector E."
	
Gibbs2Euler312::usage = 
	"E = Gibbs2Euler312[Q] translates the Gibbs
	 vector Q into the (3-1-2) Euler angle vector E."
	
Gibbs2Euler313::usage = 
	"E = Gibbs2Euler313[Q] translates the Gibbs
	 vector Q into the (3-1-3) Euler angle vector E."
	
Gibbs2Euler212::usage = 
	"E = Gibbs2Euler212[Q] translates the Gibbs
	 vector Q into the (2-1-2) Euler angle vector E."
	
Gibbs2Euler213::usage = 
	"E = Gibbs2Euler213[Q] translates the Gibbs
	 vector Q into the (2-1-3) Euler angle vector E."
	
Gibbs2Euler232::usage = 
	"E = Gibbs2Euler232[Q] translates the Gibbs
	 vector Q into the (2-3-2) Euler angle vector E."
	
Gibbs2Euler231::usage = 
	"E = Gibbs2Euler231[Q] translates the Gibbs
	 vector Q into the (2-3-1) Euler angle vector E."

Gibbs2MRP::usage =
	"Q = Gibbs2MRP[Q1] translates the Gibbs vector Q1
	into the MRP vector Q."

Gibbs2PRV::usage = 
	"Q = Gibbs2PRV[Q1] translates the Gibbs vector Q1
	into the principal rotation vector Q."

MRP2C::usage =
	"C = MRP2C[Q] returns the direction cosine 
	matrix in terms of the 3x1 MRP vector Q."

MRP2EP::usage = 
	"Q = MRP2EP[Q1] translates the MRP vector Q1
	into the Euler parameter vector Q."

MRP2Euler121::usage =
	"E = MRP2Euler121[Q] translates the MRP
	 vector Q into the (1-2-1) Euler angle vector E."

MRP2Euler123::usage =
	"E = MRP2Euler123[Q] translates the MRP
	 vector Q into the (1-2-3) Euler angle vector E."

MRP2Euler131::usage =
	"E = MRP2Euler131[Q] translates the MRP
	 vector Q into the (1-3-1) Euler angle vector E."

MRP2Euler132::usage =
	"E = MRP2Euler132[Q] translates the MRP
	 vector Q into the (1-3-2) Euler angle vector E."

MRP2Euler232::usage =
	"E = MRP2Euler232[Q] translates the MRP
	 vector Q into the (2-3-2) Euler angle vector E."

MRP2Euler231::usage =
	"E = MRP2Euler231[Q] translates the MRP
	 vector Q into the (2-3-1) Euler angle vector E."

MRP2Euler212::usage =
	"E = MRP2Euler212[Q] translates the MRP
	 vector Q into the (2-1-2) Euler angle vector E."

MRP2Euler213::usage =
	"E = MRP2Euler213[Q] translates the MRP
	 vector Q into the (2-1-3) Euler angle vector E."

MRP2Euler313::usage =
	"E = MRP2Euler313[Q] translates the MRP
	 vector Q into the (3-1-3) Euler angle vector E."

MRP2Euler323::usage =
	"E = MRP2Euler323[Q] translates the MRP
	 vector Q into the (3-2-3) Euler angle vector E."

MRP2Euler312::usage =
	"E = MRP2Euler312[Q] translates the MRP
	 vector Q into the (3-1-2) Euler angle vector E."

MRP2Euler321::usage =
	"E = MRP2Euler321[Q] translates the MRP
	 vector Q into the (3-2-1) Euler angle vector E."

MRP2Gibbs::usage = 
	"Q = MRP2Gibbs[Q1] translates the MRP vector Q1
	into the Gibbs vector Q."

MRP2PRV::usage = 
	"Q = MRP2PRV[Q1] translates the MRP vector Q1
	into the principal rotation vector Q."

MRPswitch::usage = 
	"S = MRPswitch[Q,s2] checks to see if Norm(Q) 
	is larger than s2.  If yes, then the MRP 
	vector Q is mapped to its shadow set."

PRV2C::usage = 
	"C = PRV2C[Q] returns the direction cosine 
	matrix in terms of the 3x1 principal rotation 
	vector Q. "

PRV2elem::usage = 
	"Q = PRV2elem[R] translates a prinicpal rotation 
	vector R into the corresponding principal 
	rotation element set Q."

PRV2EP::usage =
	"Q = PRV2EP[Q1] translates the principal 
	rotation vector Q into the Euler parameter 
	vector Q."

PRV2Euler121::usage = 
	"E = PRV2Euler121[Q] translates the principal 
	rotation vector Q into the (1-2-1) Euler 
	angle vector E."

PRV2Euler123::usage = 
	"E = PRV2Euler123[Q] translates the principal 
	rotation vector Q into the (1-2-3) Euler 
	angle vector E."

PRV2Euler131::usage = 
	"E = PRV2Euler131[Q] translates the principal 
	rotation vector Q into the (1-3-1) Euler 
	angle vector E."

PRV2Euler132::usage = 
	"E = PRV2Euler132[Q] translates the principal 
	rotation vector Q into the (1-3-2) Euler 
	angle vector E."

PRV2Euler232::usage = 
	"E = PRV2Euler232[Q] translates the principal 
	rotation vector Q into the (2-3-2) Euler 
	angle vector E."

PRV2Euler231::usage = 
	"E = PRV2Euler231[Q] translates the principal 
	rotation vector Q into the (2-3-1) Euler 
	angle vector E."

PRV2Euler212::usage = 
	"E = PRV2Euler212[Q] translates the principal 
	rotation vector Q into the (2-1-2) Euler 
	angle vector E."

PRV2Euler213::usage = 
	"E = PRV2Euler213[Q] translates the principal 
	rotation vector Q into the (2-1-3) Euler 
	angle vector E."

PRV2Euler313::usage = 
	"E = PRV2Euler313[Q] translates the principal 
	rotation vector Q into the (3-1-3) Euler 
	angle vector E."

PRV2Euler312::usage = 
	"E = PRV2Euler312[Q] translates the principal 
	rotation vector Q into the (3-1-2) Euler 
	angle vector E."

PRV2Euler323::usage = 
	"E = PRV2Euler323[Q] translates the principal 
	rotation vector Q into the (3-2-3) Euler 
	angle vector E."

PRV2Euler321::usage = 
	"E = PRV2Euler321[Q] translates the principal 
	rotation vector Q into the (3-2-1) Euler 
	angle vector E."

PRV2Gibbs::usage = 
	"Q = PRV2Gibbs[Q1] translates the principal 
	rotation vector Q1 into the Gibbs vector Q."

PRV2MRP::usage = 
	"Q = PRV2MRP[Q1] translates the principal 
	rotation vector Q1 into the MRP vector Q."

subEP::usage = 
	"Q = subEP[B1,B2] provides the Euler 
	parameter vector which corresponds to 
	relative rotation from B2 to B1."

subEuler121::usage = 
	"E2 = subEuler121[E,E1] computes the relative
	(1-2-1) Euler angle vector from E1 to E."

subEuler131::usage = 
	"E2 = subEuler131[E,E1] computes the relative
	(1-3-1) Euler angle vector from E1 to E."

subEuler212::usage = 
	"E2 = subEuler212[E,E1] computes the relative
	(2-1-2) Euler angle vector from E1 to E."

subEuler232::usage = 
	"E2 = subEuler232[E,E1] computes the relative
	(2-3-2) Euler angle vector from E1 to E."

subEuler313::usage = 
	"E2 = subEuler313[E,E1] computes the relative
	(3-1-3) Euler angle vector from E1 to E."

subEuler323::usage = 
	"E2 = subEuler323[E,E1] computes the relative
	(3-2-3) Euler angle vector from E1 to E."

subEuler123::usage = 
	"E2 = subEuler123[E,E1] computes the relative
	(1-2-3) Euler angle vector from E1 to E."

subEuler132::usage = 
	"E2 = subEuler132[E,E1] computes the relative
	(1-3-2) Euler angle vector from E1 to E."

subEuler213::usage = 
	"E2 = subEuler213[E,E1] computes the relative
	(2-1-3) Euler angle vector from E1 to E."

subEuler231::usage = 
	"E2 = subEuler231[E,E1] computes the relative
	(2-3-1) Euler angle vector from E1 to E."

subEuler312::usage = 
	"E2 = subEuler312[E,E1] computes the relative
	(3-1-2) Euler angle vector from E1 to E."

subEuler321::usage = 
	"E2 = subEuler321[E,E1] computes the relative
	(3-2-1) Euler angle vector from E1 to E."

subGibbs::usage = 
	"Q = subGibbs[Q1,Q2] provides the Gibbs vector
	which corresponds to relative rotation from Q2
	to Q1."

subMRP::usage = 
	"Q = subMRP[Q1,Q2] provides the MRP vector
	which corresponds to relative rotation from Q2
	to Q1."

subPRV::usage = 
	"Q = subPRV(Q1,Q2) provides the prinipal 
	rotation vector which corresponds to 
	relative principal rotation from Q2 to Q1."

tilde::usage = 
	"W = tilde(w) returns the 3x3 tilde matrix of the 3x1 vector w."

























Picheck::usage = 
	"Makes sure that the angle lies within +/- Pi."

Begin["`Private`"]

addEP[b1_,b2_]:={
	b2[[1]]*b1[[1]]-b2[[2]]*b1[[2]]-b2[[3]]*b1[[3]]-b2[[4]]*b1[[4]],
	b2[[2]]*b1[[1]]+b2[[1]]*b1[[2]]+b2[[4]]*b1[[3]]-b2[[3]]*b1[[4]],
	b2[[3]]*b1[[1]]-b2[[4]]*b1[[2]]+b2[[1]]*b1[[3]]+b2[[2]]*b1[[4]],
	b2[[4]]*b1[[1]]+b2[[3]]*b1[[2]]-b2[[2]]*b1[[3]]+b2[[1]]*b1[[4]]}


addEuler121[e1_,e2_]:= Block[{cp1,cp2,sp1,sp2,dum},
	cp1 = Cos[e1[[2]]];
	cp2 = Cos[e2[[2]]];
	sp1 = Sin[e1[[2]]];
	sp2 = Sin[e2[[2]]];
	dum = e1[[3]]+e2[[1]];
	
	q = {1,1,1};
	q[[2]] = ArcCos[cp1*cp2-sp1*sp2*Cos[dum]];
	cp3 = Cos[q[[2]]];
	q[[1]] = Picheck[e1[[1]] + ArcTan[cp2-cp3*cp1,sp1*sp2*Sin[dum]]];
	q[[3]] = Picheck[e2[[3]] + ArcTan[cp1-cp3*cp2,sp1*sp2*Sin[dum]]]; 
	q
	]


addEuler123[e1_,e2_]:=Block[{C1,C2,C3},
	C1 = Euler1232C[e1];
	C2 = Euler1232C[e2];
	C3 = C2.C1;
	C2Euler123[C3]
	]

addEuler131[e1_,e2_]:= Block[{cp1,cp2,sp1,sp2,dum},
	cp1 = Cos[e1[[2]]];
	cp2 = Cos[e2[[2]]];
	sp1 = Sin[e1[[2]]];
	sp2 = Sin[e2[[2]]];
	dum = e1[[3]]+e2[[1]];
	
	q = {1,1,1};
	q[[2]] = ArcCos[cp1*cp2-sp1*sp2*Cos[dum]];
	cp3 = Cos[q[[2]]];
	q[[1]] = Picheck[e1[[1]] + ArcTan[cp2-cp3*cp1,sp1*sp2*Sin[dum]]];
	q[[3]] = Picheck[e2[[3]] + ArcTan[cp1-cp3*cp2,sp1*sp2*Sin[dum]]]; 
	q
	]

addEuler132[e1_,e2_]:=Block[{C1,C2,Cm},
	C1 = Euler1322C[e1];
	C2 = Euler1322C[e2];
	Cm = C2.C1;
	C2Euler132[Cm]
	]

addEuler212[e1_,e2_]:= Block[{cp1,cp2,sp1,sp2,dum},
	cp1 = Cos[e1[[2]]];
	cp2 = Cos[e2[[2]]];
	sp1 = Sin[e1[[2]]];
	sp2 = Sin[e2[[2]]];
	dum = e1[[3]]+e2[[1]];
	
	q = {1,1,1};
	q[[2]] = ArcCos[cp1*cp2-sp1*sp2*Cos[dum]];
	cp3 = Cos[q[[2]]];
	q[[1]] = Picheck[e1[[1]] + ArcTan[cp2-cp3*cp1,sp1*sp2*Sin[dum]]];
	q[[3]] = Picheck[e2[[3]] + ArcTan[cp1-cp3*cp2,sp1*sp2*Sin[dum]]]; 
	q
	]

addEuler213[e1_,e2_]:=Block[{C1,C2,Cm},
	C1 = Euler2132C[e1];
	C2 = Euler2132C[e2];
	Cm = C2.C1;
	C2Euler213[Cm]
	]

addEuler231[e1_,e2_]:=Block[{C1,C2,Cm},
	C1 = Euler2312C[e1];
	C2 = Euler2312C[e2];
	Cm = C2.C1;
	C2Euler231[Cm]
	]

addEuler232[e1_,e2_]:= Block[{cp1,cp2,sp1,sp2,dum},
	cp1 = Cos[e1[[2]]];
	cp2 = Cos[e2[[2]]];
	sp1 = Sin[e1[[2]]];
	sp2 = Sin[e2[[2]]];
	dum = e1[[3]]+e2[[1]];
	
	q = {1,1,1};
	q[[2]] = ArcCos[cp1*cp2-sp1*sp2*Cos[dum]];
	cp3 = Cos[q[[2]]];
	q[[1]] = Picheck[e1[[1]] + ArcTan[cp2-cp3*cp1,sp1*sp2*Sin[dum]]];
	q[[3]] = Picheck[e2[[3]] + ArcTan[cp1-cp3*cp2,sp1*sp2*Sin[dum]]]; 
	q
	]

addEuler312[e1_,e2_]:=Block[{C1,C2,Cm},
	C1 = Euler3122C[e1];
	C2 = Euler3122C[e2];
	Cm = C2.C1;
	C2Euler312[Cm]
	]

addEuler313[e1_,e2_]:= Block[{cp1,cp2,sp1,sp2,dum},
	cp1 = Cos[e1[[2]]];
	cp2 = Cos[e2[[2]]];
	sp1 = Sin[e1[[2]]];
	sp2 = Sin[e2[[2]]];
	dum = e1[[3]]+e2[[1]];
	
	q = {1,1,1};
	q[[2]] = ArcCos[cp1*cp2-sp1*sp2*Cos[dum]];
	cp3 = Cos[q[[2]]];
	q[[1]] = Picheck[e1[[1]] + ArcTan[cp2-cp3*cp1,sp1*sp2*Sin[dum]]];
	q[[3]] = Picheck[e2[[3]] + ArcTan[cp1-cp3*cp2,sp1*sp2*Sin[dum]]]; 
	q
	]

addEuler321[e1_,e2_]:=Block[{C1,C2,Cm},
	C1 = Euler3212C[e1];
	C2 = Euler3212C[e2];
	Cm = C2.C1;
	C2Euler321[Cm]
	]

addEuler323[e1_,e2_]:= Block[{cp1,cp2,sp1,sp2,dum},
	cp1 = Cos[e1[[2]]];
	cp2 = Cos[e2[[2]]];
	sp1 = Sin[e1[[2]]];
	sp2 = Sin[e2[[2]]];
	dum = e1[[3]]+e2[[1]];
	
	q = {1,1,1};
	q[[2]] = ArcCos[cp1*cp2-sp1*sp2*Cos[dum]];
	cp3 = Cos[q[[2]]];
	q[[1]] = Picheck[e1[[1]] + ArcTan[cp2-cp3*cp1,sp1*sp2*Sin[dum]]];
	q[[3]] = Picheck[e2[[3]] + ArcTan[cp1-cp3*cp2,sp1*sp2*Sin[dum]]]; 
	q
	]

addGibbs[q1_,q2_] := (q1+q2+Cross[q1,q2])/(1-Dot[q1,q2])

addMRP[q1_,q2_] := ((1-q1.q1)*q2+(1-q2.q2)*q1+
	2*Cross[q1,q2])/(1+q1.q1*q2.q2-2*q1.q2)

addPRV[q1v_,q2v_] := Block[{q1,q2,cp1,cp2,sp1,sp2,e1,e2,p,sp,e},
	q1 = PRV2elem[q1v];
	q2 = PRV2elem[q2v];
	cp1 = Cos[q1[[1]]/2];
	cp2 = Cos[q2[[1]]/2]; 
	sp1 = Sin[q1[[1]]/2];
	sp2 = Sin[q2[[1]]/2];
	e1 = {q1[[2]],q1[[3]],q1[[4]]};
	e2 = {q2[[2]],q2[[3]],q2[[4]]};
	
	p = 2*ArcCos[cp1*cp2-sp1*sp2*e1.e2];
	sp = Sin[p/2];
	e = (cp1*sp2*e2+cp2*sp1*e1+sp1*sp2*Cross[e1,e2])/sp;
	p*e
	]

BinvEP[q_]:={{-q[[2]],q[[1]],q[[4]],-q[[3]]},
	{-q[[3]],-q[[4]],q[[1]],q[[2]]},
	{-q[[4]],q[[3]],-q[[2]],q[[1]]}}

BinvEuler121[q_]:=Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{c2,0,1},{s2*s3,c3,0},{s2*c3,-s3,0}}
	]

BinvEuler123[q_]:=Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{c2*c3,s3,0},{-c2*s3,c3,0},{s2,0,1}}
	]

BinvEuler131[q_]:=Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{c2,0,1},{-s2*c3,s3,0},{s2*s3,c3,0}}
	]

BinvEuler132[q_]:=Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{c2*c3,-s3,0},{-s2,0,1},{c2*s3,c3,0}}
	]

BinvEuler212[q_]:=Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{s2*s3,c3,0},{c2,0,1},{-s2*c3,s3,0}}
	]

BinvEuler213[q_]:=Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{c2*s3,c3,0},{c2*c3,-s3,0},{-s2,0,1}}
	]

BinvEuler231[q_]:=Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{s2,0,1},{c2*c3,s3,0},{-c2*s3,c3,0}}
	]

BinvEuler232[q_]:=Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{s2*c3,-s3,0},{c2,0,1},{s2*s3,c3,0}}
	]

BinvEuler312[q_]:=Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{-c2*s3,c3,0},{s2,0,1},{c2*c3,s3,0}}
	]

BinvEuler313[q_]:=Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{s2*s3,c3,0},{s2*c3,-s3,0},{c2,0,1}}
	]

BinvEuler321[q_]:=Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{-s2,0,1},{c2*s3,c3,0},{c2*c3,-s3,0}}
	]

BinvEuler323[q_]:=Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{-s2*c3,s3,0},{s2*s3,c3,0},{c2,0,1}}
	]

BinvGibbs[q_]:={{1,q[[3]],-q[[2]]},
{-q[[3]],1,q[[1]]},{q[[2]],-q[[1]],1}}/(1+q.q)
	
BinvMRP[q_] := Block[{s2},
	s2 = q.q;
	{{1-s2+2*q[[1]]*q[[1]],2*(q[[1]]*q[[2]]+q[[3]]),
	2*(q[[1]]*q[[3]]-q[[2]])},
	{2*(q[[2]]*q[[1]]-q[[3]]),1-s2+2*q[[2]]*q[[2]],
	2*(q[[2]]*q[[3]]+q[[1]])},
	{2*(q[[3]]*q[[1]]+q[[2]]),2*(q[[3]]*q[[2]]-q[[1]]),
	1-s2+2*q[[3]]*q[[3]]}}/(1+s2)/(1+s2)
	]

BinvPRV[q_] := Block[{p,c1,c2},
	p = Sqrt[q.q];
	c1 = (1-Cos[p])/p/p;
	c2 = (p-Sin[p])/p/p/p;
	
	{{1-c2*(q[[2]]*q[[2]]+q[[3]]*q[[3]]),
	c1*q[[3]]+c2*q[[1]]*q[[2]],
	-c1*q[[2]]+c2*q[[1]]*q[[3]]},
	{-c1*q[[3]]+c2*q[[1]]*q[[2]],
	1-c2*(q[[1]]*q[[1]]+q[[3]]*q[[3]]),
	c1*q[[1]]+c2*q[[2]]*q[[3]]},
	{c1*q[[2]]+c2*q[[3]]*q[[1]],
	-c1*q[[1]]+c2*q[[3]]*q[[2]],
	1-c2*(q[[1]]*q[[1]]+q[[2]]*q[[2]])}}
	]

BmatEP[q_] := {{-q[[2]],-q[[3]],-q[[4]]},
{q[[1]],-q[[4]],q[[3]]},
{q[[4]],q[[1]],-q[[2]]},
{-q[[3]],q[[2]],q[[1]]}}
	
BmatEuler121[q_] := Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{0,s3,c3},{0,s2*c3,-s2*s3},
	{s2,-c2*s3,-c2*c3}}/s2
	]

BmatEuler123[q_] := Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{c3,-s3,0},{c2*s3,c2*c3,0},
	{-s2*c3,s2*s3,c2}}/c2
	]

BmatEuler131[q_] := Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{0,-c3,s3},{0,s2*s3,s2*c3},{s2,c2*c3,-c2*s3}}/s2
	]

BmatEuler132[q_] := Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{c3,0,s3},{-c2*s3,0,c2*c3},{s2*c3,c2,s2*s3}}/c2
	]

BmatEuler212[q_] := Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{s3,0,-c3},{s2*c3,0,s2*s3},{-c2*s3,s2,c2*c3}}/s2
	]

BmatEuler213[q_] := Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{s3,c3,0},{c2*c3,-c2*s3,0},{s2*s3,s2*c3,c2}}/c2
	]

BmatEuler231[q_] := Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{0,c3,-s3},{0,c2*s3,c2*c3},{c2,-s2*c3,s2*s3}}/c2
	]

BmatEuler232[q_] := Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{c3,0,s3},{-s2*s3,0,s2*c3},{-c2*c3,s2,-c2*s3}}/s2
	]

BmatEuler312[q_] := Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{-s3,0,c3},{c2*c3,0,c2*s3},{s2*s3,c2,-s2*c3}}/c2
	]

BmatEuler313[q_] := Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{s3,c3,0},{c3*s2,-s3*s2,0},{-s3*c2,-c3*c2,s2}}/s2
	]

BmatEuler321[q_] := Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{0,s3,c3},{0,c2*c3,-c2*s3},{c2,s2*s3,s2*c3}}/c2
	]

BmatEuler323[q_] := Block[{s2,c2,s3,c3},
	s2 = Sin[q[[2]]];
	c2 = Cos[q[[2]]];
	s3 = Sin[q[[3]]];
	c3 = Cos[q[[3]]];
	
	{{-c3,s3,0},{s2*s3,s2*c3,0},{c2*c3,-c2*s3,s2}}/s2
	]

BmatGibbs[q_] := {{1+q[[1]]*q[[1]],
q[[1]]*q[[2]]-q[[3]],q[[1]]*q[[3]]+q[[2]]},
{q[[2]]*q[[1]]+q[[3]],1+q[[2]]*q[[2]],
q[[2]]*q[[3]]-q[[1]]},
{q[[3]]*q[[1]]-q[[2]],q[[3]]*q[[2]]+q[[1]],
1+q[[3]]*q[[3]]}}

BmatMRP[q_] := Block[{s2},
	s2 = q.q;
	{{1-s2+2*q[[1]]*q[[1]],2*(q[[1]]*q[[2]]-q[[3]]),
	2*(q[[1]]*q[[3]]+q[[2]])},
	{2*(q[[2]]*q[[1]]+q[[3]]),1-s2+2*q[[2]]*q[[2]],
	2*(q[[2]]*q[[3]]-q[[1]])},
	{2*(q[[3]]*q[[1]]-q[[2]]),2*(q[[3]]*q[[2]]+q[[1]]),
	1-s2+2*q[[3]]*q[[3]]}}
	]

BmatPRV[q_] := Block[{p,c},
	p = Sqrt[q.q];
	c = 1/p/p*(1-p/2*Cot[p/2]);
	
	{{1-c*(q[[2]]*q[[2]]+q[[3]]*q[[3]]),
	-q[[3]]/2+c*q[[1]]*q[[2]],
	q[[2]]/2+c*q[[1]]*q[[3]]},
	{q[[3]]/2+c*q[[1]]*q[[2]],
	1-c*(q[[1]]*q[[1]]+q[[3]]*q[[3]]),
	-q[[1]]/2+c*q[[2]]*q[[3]]},
	{-q[[2]]/2+c*q[[1]]*q[[3]],
	q[[1]]/2+c*q[[2]]*q[[3]],
	1-c*(q[[1]]*q[[1]]+q[[2]]*q[[2]])}}
	]

C2EP[Cm_] := Block[{tr,b,b2,bi,largest},
	tr = Cm[[1,1]]+Cm[[2,2]]+Cm[[3,3]];
	b2 = {(1+tr)/4,
	(1+2*Cm[[1,1]]-tr)/4,
	(1+2*Cm[[2,2]]-tr)/4,
	(1+2*Cm[[3,3]]-tr)/4};
	
	largest = Max[b2];
	If[largest==b2[[1]],{
		bi = Sqrt[b2[[1]]];
		b = {bi,
			(Cm[[2,3]]-Cm[[3,2]])/4/bi,
			(Cm[[3,1]]-Cm[[1,3]])/4/bi,
			(Cm[[1,2]]-Cm[[2,1]])/4/bi};
		},
		If[largest == b2[[2]],{
			bi = Sqrt[b2[[2]]];
			b = {(Cm[[2,3]]-Cm[[3,2]])/4/bi,
				bi,
				(Cm[[1,2]]+Cm[[2,1]])/4/bi,
				(Cm[[3,1]]+Cm[[1,3]])/4/bi};
				If[b[[1]]<0,b=-b];
			},
			If[largest == b2[[3]],{
				bi = Sqrt[b2[[3]]];
				b = {(Cm[[3,1]]-Cm[[1,3]])/4/bi,
					(Cm[[1,2]]+Cm[[2,1]])/4/bi,
					bi,
					(Cm[[2,3]]+Cm[[3,2]])/4/bi};
				If[b[[1]]<0,b=-b];
				},
				If[largest == b2[[4]],{
					bi = Sqrt[b2[[4]]];
					b = {(Cm[[1,2]]-Cm[[2,1]])/4/bi,
					(Cm[[3,1]]+Cm[[1,3]])/4/bi,
					(Cm[[2,3]]+Cm[[3,2]])/4/bi,bi};
					If[b[[1]]<0,b=-b];
					}];
				];
			];
		];
	b
	]

C2Euler121[Cm_]:={ArcTan[-Cm[[1,3]],Cm[[1,2]]],
	ArcCos[Cm[[1,1]]],
	ArcTan[Cm[[3,1]],Cm[[2,1]]]}

C2Euler123[Cm_] := {ArcTan[Cm[[3,3]],-Cm[[3,2]]],
	ArcSin[Cm[[3,1]]],
	ArcTan[Cm[[1,1]],-Cm[[2,1]]]}
	
C2Euler131[Cm_] := {ArcTan[Cm[[1,2]],Cm[[1,3]]],
	ArcCos[Cm[[1,1]]],
	ArcTan[-Cm[[2,1]],Cm[[3,1]]]}

C2Euler132[Cm_] := {ArcTan[Cm[[2,2]],Cm[[2,3]]],
	ArcSin[-Cm[[2,1]]],
	ArcTan[Cm[[1,1]],Cm[[3,1]]]}

C2Euler212[Cm_] := {ArcTan[Cm[[2,3]],Cm[[2,1]]],
	ArcCos[Cm[[2,2]]],
	ArcTan[-Cm[[3,2]],Cm[[1,2]]]}

C2Euler213[Cm_] := {ArcTan[Cm[[3,3]],Cm[[3,1]]],
	ArcSin[-Cm[[3,2]]],
	ArcTan[Cm[[2,2]],Cm[[1,2]]]}

C2Euler231[Cm_] := {ArcTan[Cm[[1,1]],-Cm[[1,3]]],
	ArcSin[Cm[[1,2]]],
	ArcTan[Cm[[2,2]],-Cm[[3,2]]]}

C2Euler232[Cm_] := {ArcTan[-Cm[[2,1]],Cm[[2,3]]],
	ArcCos[Cm[[2,2]]],
	ArcTan[Cm[[1,2]],Cm[[3,2]]]}

C2Euler312[Cm_] := {ArcTan[Cm[[2,2]],-Cm[[2,1]]],
	ArcSin[Cm[[2,3]]],
	ArcTan[Cm[[3,3]],-Cm[[1,3]]]}

C2Euler313[Cm_] := {ArcTan[-Cm[[3,2]],Cm[[3,1]]],
	ArcCos[Cm[[3,3]]],
	ArcTan[Cm[[2,3]],Cm[[1,3]]]}

C2Euler321[Cm_] := {ArcTan[Cm[[1,1]],Cm[[1,2]]],
	ArcSin[-Cm[[1,3]]],
	ArcTan[Cm[[3,3]],Cm[[2,3]]]}

C2Euler323[Cm_] := {ArcTan[Cm[[3,1]],Cm[[3,2]]],
	ArcCos[Cm[[3,3]]],
	ArcTan[-Cm[[1,3]],Cm[[2,3]]]}

C2Gibbs[Cm_] := Block[{b},
	b = C2EP[Cm];
	{b[[2]],b[[3]],b[[4]]}/b[[1]]
	]

C2MRP[Cm_] := Block[{b},
	b = C2EP[Cm];
	{b[[2]],b[[3]],b[[4]]}/(1+b[[1]])
	]

C2PRV[Cm_] := Block[{cp,p,sp},
	cp = (Cm[[1,1]]+Cm[[2,2]]+Cm[[3,3]]-1)/2;
	p = ArcCos[cp];
	sp = p/2/Sin[p];
	{Cm[[2,3]]-Cm[[3,2]],
	Cm[[3,1]]-Cm[[1,3]],
	Cm[[1,2]]-Cm[[2,1]]}*sp
	]

dEP[q_,w_] := 1/2*BmatEP[q].w;

dEuler121[q_,w_] := BmatEuler121[q].w;

dEuler123[q_,w_] := BmatEuler123[q].w;

dEuler131[q_,w_] := BmatEuler131[q].w;

dEuler132[q_,w_] := BmatEuler132[q].w;

dEuler212[q_,w_] := BmatEuler212[q].w;

dEuler213[q_,w_] := BmatEuler213[q].w;

dEuler231[q_,w_] := BmatEuler231[q].w;

dEuler232[q_,w_] := BmatEuler232[q].w;

dEuler312[q_,w_] := BmatEuler312[q].w;

dEuler313[q_,w_] := BmatEuler313[q].w;

dEuler321[q_,w_] := BmatEuler321[q].w;

dEuler323[q_,w_] := BmatEuler323[q].w;

dGibbs[q_,w_] := 1/2*BmatGibbs[q].w;

dMRP[q_,w_] := 1/4*BmatMRP[q].w;

dPRV[q_,w_] := BmatPRV[q].w;

elem2PRV[r_] := {r[[2]]*r[[1]],r[[3]]*r[[1]],r[[4]]*r[[1]]};

EP2C[q_] := Block[{q0,q1,q2,q3,Cm},
	q0 = q[[1]];
	q1 = q[[2]];
	q2 = q[[3]];
	q3 = q[[4]];
	
	Cm = IdentityMatrix[3];
	Cm[[1,1]] = q0*q0+q1*q1-q2*q2-q3*q3;
	Cm[[1,2]] = 2*(q1*q2+q0*q3);
	Cm[[1,3]] = 2*(q1*q3-q0*q2);
	Cm[[2,1]] = 2*(q1*q2-q0*q3);
	Cm[[2,2]] = q0*q0-q1*q1+q2*q2-q3*q3;
	Cm[[2,3]] = 2*(q2*q3+q0*q1);
	Cm[[3,1]] = 2*(q1*q3 + q0*q2);
	Cm[[3,2]] = 2*(q2*q3-q0*q1);
	Cm[[3,3]] = q0*q0-q1*q1-q2*q2+q3*q3;
	Cm
	]

EP2Euler121[q_] := Block[{t1,t2,e},
	t1 = ArcTan[q[[3]],q[[4]]];
	t2 = ArcTan[q[[1]],q[[2]]];
	
	e = {t1+t2,
		2*ArcCos[Sqrt[q[[1]]*q[[1]]+q[[2]]*q[[2]]]],
		t2-t1}
	]

EP2Euler123[q_] := Block[{q0,q1,q2,q3,e},
	q0 = q[[1]];
	q1 = q[[2]];
	q2 = q[[3]];
	q3 = q[[4]];
	
	e = {ArcTan[q0*q0-q1*q1-q2*q2+q3*q3,-2*(q2*q3-q0*q1)],
		ArcSin[2*(q1*q3 + q0*q2)],
		ArcTan[q0*q0+q1*q1-q2*q2-q3*q3,-2*(q1*q2-q0*q3)]}
	]

EP2Euler131[q_] := Block[{t1,t2,e},
	t1 = ArcTan[q[[4]],q[[3]]];
	t2 = ArcTan[q[[1]],q[[2]]];
	
	e = {t2-t1,
		2*ArcCos[Sqrt[q[[1]]*q[[1]]+q[[2]]*q[[2]]]],
		t2+t1}
	]

EP2Euler132[q_] := Block[{q0,q1,q2,q3,e},
	q0 = q[[1]];
	q1 = q[[2]];
	q2 = q[[3]];
	q3 = q[[4]];
	
	e = {ArcTan[q0*q0-q1*q1+q2*q2-q3*q3,2*(q2*q3+q0*q1)],
		ArcSin[-2*(q1*q2-q0*q3)],
		ArcTan[q0*q0+q1*q1-q2*q2-q3*q3,2*(q1*q3 + q0*q2)]}
	]

EP2Euler212[q_] := Block[{t1,t2,e},
	t1 = ArcTan[q[[2]],q[[4]]];
	t2 = ArcTan[q[[1]],q[[3]]];
	
	e = {t2-t1,
		2*ArcCos[Sqrt[q[[1]]*q[[1]]+q[[3]]*q[[3]]]],
		t2+t1}
	]

EP2Euler213[q_] := Block[{q0,q1,q2,q3,e},
	q0 = q[[1]];
	q1 = q[[2]];
	q2 = q[[3]];
	q3 = q[[4]];
	
	e = {ArcTan[q0*q0-q1*q1-q2*q2+q3*q3,2*(q1*q3 + q0*q2)],
		ArcSin[-2*(q2*q3-q0*q1)],
		ArcTan[q0*q0-q1*q1+q2*q2-q3*q3,2*(q1*q2+q0*q3)]}
	]

EP2Euler231[q_] := Block[{q0,q1,q2,q3,e},
	q0 = q[[1]];
	q1 = q[[2]];
	q2 = q[[3]];
	q3 = q[[4]];
	
	e = {ArcTan[ q0*q0+q1*q1-q2*q2-q3*q3,-2*(q1*q3-q0*q2)],
		ArcSin[2*(q1*q2+q0*q3)],
		ArcTan[q0*q0-q1*q1+q2*q2-q3*q3,-2*(q2*q3-q0*q1)]}
	]

EP2Euler232[q_] := Block[{t1,t2,e},
	t1 = ArcTan[q[[4]],q[[2]]];
	t2 = ArcTan[q[[1]],q[[3]]];
	
	e = {t2+t1,
		2*ArcCos[Sqrt[q[[1]]*q[[1]]+q[[3]]*q[[3]]]],
		t2-t1}
	]

EP2Euler312[q_] := Block[{q0,q1,q2,q3,e},
	q0 = q[[1]];
	q1 = q[[2]];
	q2 = q[[3]];
	q3 = q[[4]];
	
	e = {ArcTan[q0*q0-q1*q1+q2*q2-q3*q3,-2*(q1*q2-q0*q3)],
		ArcSin[2*(q2*q3+q0*q1)],
		ArcTan[q0*q0-q1*q1-q2*q2+q3*q3,-2*(q1*q3-q0*q2)]}
	]

EP2Euler313[q_] := Block[{t1,t2,e},
	t1 = ArcTan[q[[2]],q[[3]]];
	t2 = ArcTan[q[[1]],q[[4]]];
	
	e = {t2+t1,
		2*ArcCos[Sqrt[q[[1]]*q[[1]]+q[[4]]*q[[4]]]],
		t2-t1}
	]

EP2Euler321[q_] := Block[{q0,q1,q2,q3,e},
	q0 = q[[1]];
	q1 = q[[2]];
	q2 = q[[3]];
	q3 = q[[4]];
	
	e = {ArcTan[q0*q0+q1*q1-q2*q2-q3*q3,2*(q1*q2+q0*q3)],
		ArcSin[-2*(q1*q3-q0*q2)],
		ArcTan[q0*q0-q1*q1-q2*q2+q3*q3,2*(q2*q3+q0*q1)]}
	]

EP2Euler323[q_] := Block[{t1,t2,e},
	t1 = ArcTan[q[[3]],q[[2]]];
	t2 = ArcTan[q[[1]],q[[4]]];
	
	e = {t2-t1,
		2*ArcCos[Sqrt[q[[1]]*q[[1]]+q[[4]]*q[[4]]]],
		t2+t1}
	]

EP2Gibbs[q_] := {q[[2]]/q[[1]],
				q[[3]]/q[[1]],
				q[[4]]/q[[1]]};

EP2MRP[q_] := {q[[2]]/(1+q[[1]]),
				q[[3]]/(1+q[[1]]),
				q[[4]]/(1+q[[1]])};

EP2PRV[q_] := Block[{p,sp},
	p = 2*ArcCos[q[[1]]];
	sp = Sin[p/2];
	{q[[2]]/sp*p,q[[3]]/sp*p,q[[4]]/sp*p}
	]

Euler1[x_] := {{1,0,0},{0,Cos[x],Sin[x]},{0,-Sin[x],Cos[x]}}

Euler2[x_] := {{Cos[x],0,-Sin[x]},{0,1,0},{Sin[x],0,Cos[x]}}

Euler3[x_] := {{Cos[x],Sin[x],0},{-Sin[x],Cos[x],0},{0,0,1}}

Euler1212C[q_] := Block[{st1,ct1,st2,ct2,st3,ct3},
	st1 = Sin[q[[1]]];
	ct1 = Cos[q[[1]]];
	st2 = Sin[q[[2]]];
	ct2 = Cos[q[[2]]];
	st3 = Sin[q[[3]]];
	ct3 = Cos[q[[3]]];
	
	{{ct2,st1*st2,-ct1*st2},
	{st2*st3,ct1*ct3-ct2*st1*st3,ct3*st1+ct1*ct2*st3},
	{ct3*st2,-ct2*ct3*st1-ct1*st3,ct1*ct2*ct3-st1*st3}}
	]
	
Euler1232C[q_] := Block[{st1,ct1,st2,ct2,st3,ct3},
	st1 = Sin[q[[1]]];
	ct1 = Cos[q[[1]]];
	st2 = Sin[q[[2]]];
	ct2 = Cos[q[[2]]];
	st3 = Sin[q[[3]]];
	ct3 = Cos[q[[3]]];
	
	{{ct2*ct3,
	ct3*st1*st2+ct1*st3,
	st1*st3-ct1*ct3*st2},
	{-ct2*st3,
	ct1*ct3-st1*st2*st3,
	ct3*st1+ct1*st2*st3},
	{st2,
	-ct2*st1,
	ct1*ct2}}
	]
	
Euler1312C[q_] := Block[{st1,ct1,st2,ct2,st3,ct3},
	st1 = Sin[q[[1]]];
	ct1 = Cos[q[[1]]];
	st2 = Sin[q[[2]]];
	ct2 = Cos[q[[2]]];
	st3 = Sin[q[[3]]];
	ct3 = Cos[q[[3]]];
	
	{{ct2,
	ct1*st2,
	st1*st2},
	{-ct3*st2,
	ct1*ct2*ct3-st1*st3,
	ct2*ct3*st1+ct1*st3},
	{st2*st3,
	-ct3*st1-ct1*ct2*st3,
	ct1*ct3-ct2*st1*st3}}
	]
	
Euler1322C[q_] := Block[{st1,ct1,st2,ct2,st3,ct3},
	st1 = Sin[q[[1]]];
	ct1 = Cos[q[[1]]];
	st2 = Sin[q[[2]]];
	ct2 = Cos[q[[2]]];
	st3 = Sin[q[[3]]];
	ct3 = Cos[q[[3]]];
	
	{{ct2*ct3,
	ct1*ct3*st2+st1*st3,
	ct3*st1*st2-ct1*st3},
	{-st2,
	ct1*ct2,
	ct2*st1},
	{ct2*st3,
	-ct3*st1+ct1*st2*st3,
	ct1*ct3+st1*st2*st3}}
	]
	
Euler2122C[q_] := Block[{st1,ct1,st2,ct2,st3,ct3},
	st1 = Sin[q[[1]]];
	ct1 = Cos[q[[1]]];
	st2 = Sin[q[[2]]];
	ct2 = Cos[q[[2]]];
	st3 = Sin[q[[3]]];
	ct3 = Cos[q[[3]]];
	
	{{ct1*ct3-ct2*st1*st3,
	st2*st3,
	-ct3*st1-ct1*ct2*st3},
	{st1*st2,
	ct2,
	ct1*st2},
	{ct2*ct3*st1+ct1*st3,
	-ct3*st2,
	ct1*ct2*ct3-st1*st3}}
	]
	
Euler2132C[q_] := Block[{st1,ct1,st2,ct2,st3,ct3},
	st1 = Sin[q[[1]]];
	ct1 = Cos[q[[1]]];
	st2 = Sin[q[[2]]];
	ct2 = Cos[q[[2]]];
	st3 = Sin[q[[3]]];
	ct3 = Cos[q[[3]]];
	
	{{ct1*ct3+st1*st2*st3,
	ct2*st3,
	-ct3*st1+ct1*st2*st3},
	{ct3*st1*st2-ct1*st3,
	ct2*ct3,
	ct1*ct3*st2 + st1*st3},
	{ct2*st1,
	-st2,
	ct1*ct2}}
	]
	
Euler2312C[q_] := Block[{st1,ct1,st2,ct2,st3,ct3},
	st1 = Sin[q[[1]]];
	ct1 = Cos[q[[1]]];
	st2 = Sin[q[[2]]];
	ct2 = Cos[q[[2]]];
	st3 = Sin[q[[3]]];
	ct3 = Cos[q[[3]]];
	
	{{ct1*ct2,
	st2,
	-ct2*st1},
	{-ct1*ct3*st2+st1*st3,
	ct2*ct3,
	ct3*st1*st2+ct1*st3},
	{ct3*st1+ct1*st2*st3,
	-ct2*st3,
	ct1*ct3-st1*st2*st3}}
	]
	
Euler2322C[q_] := Block[{st1,ct1,st2,ct2,st3,ct3},
	st1 = Sin[q[[1]]];
	ct1 = Cos[q[[1]]];
	st2 = Sin[q[[2]]];
	ct2 = Cos[q[[2]]];
	st3 = Sin[q[[3]]];
	ct3 = Cos[q[[3]]];
	
	{{ct1*ct2*ct3-st1*st3,
	ct3*st2,
	-ct2*ct3*st1-ct1*st3},
	{-ct1*st2,
	ct2,
	st1*st2},
	{ct3*st1+ct1*ct2*st3,
	st2*st3,
	ct1*ct3-ct2*st1*st3}}
	]
	
Euler3122C[q_] := Block[{st1,ct1,st2,ct2,st3,ct3},
	st1 = Sin[q[[1]]];
	ct1 = Cos[q[[1]]];
	st2 = Sin[q[[2]]];
	ct2 = Cos[q[[2]]];
	st3 = Sin[q[[3]]];
	ct3 = Cos[q[[3]]];
	
	{{ct1*ct3-st1*st2*st3,
	ct3*st1+ct1*st2*st3,
	-ct2*st3},
	{-ct2*st1,
	ct1*ct2,
	st2},
	{ct3*st1*st2+ct1*st3,
	st1*st3-ct1*ct3*st2,
	ct2*ct3}}
	]
	
Euler3132C[q_] := Block[{st1,ct1,st2,ct2,st3,ct3},
	st1 = Sin[q[[1]]];
	ct1 = Cos[q[[1]]];
	st2 = Sin[q[[2]]];
	ct2 = Cos[q[[2]]];
	st3 = Sin[q[[3]]];
	ct3 = Cos[q[[3]]];
	
	{{ct3*ct1-st3*ct2*st1,
	ct3*st1+st3*ct2*ct1,
	st3*st2},
	{-st3*ct1-ct3*ct2*st1,
	-st3*st1+ct3*ct2*ct1,
	ct3*st2},
	{st2*st1,
	-st2*ct1,
	ct2}}
	]
	
Euler3212C[q_] := Block[{st1,ct1,st2,ct2,st3,ct3},
	st1 = Sin[q[[1]]];
	ct1 = Cos[q[[1]]];
	st2 = Sin[q[[2]]];
	ct2 = Cos[q[[2]]];
	st3 = Sin[q[[3]]];
	ct3 = Cos[q[[3]]];
	
	{{ct2*ct1,
	ct2*st1,
	-st2},
	{st3*st2*ct1-ct3*st1,
	st3*st2*st1+ct3*ct1,
	st3*ct2},
	{ct3*st2*ct1+st3*st1,
	ct3*st2*st1-st3*ct1,
	ct3*ct2}}
	]
	
Euler3232C[q_] := Block[{st1,ct1,st2,ct2,st3,ct3},
	st1 = Sin[q[[1]]];
	ct1 = Cos[q[[1]]];
	st2 = Sin[q[[2]]];
	ct2 = Cos[q[[2]]];
	st3 = Sin[q[[3]]];
	ct3 = Cos[q[[3]]];
	
	{{ct1*ct2*ct3-st1*st3,
	ct2*ct3*st1+ct1*st3,
	-ct3*st2},
	{-ct3*st1-ct1*ct2*st3,
	ct1*ct3-ct2*st1*st3,
	st2*st3},
	{ct1*st2,
	st1*st2,
	ct2}}
	]
	
Euler1212EP[e_] := Block[{e1,e2,e3},
	e1 = e[[1]]/2;
	e2 = e[[2]]/2;
	e3 = e[[3]]/2;
	
	{Cos[e2]*Cos[e1+e3],
	Cos[e2]*Sin[e1+e3],
	Sin[e2]*Cos[e1-e3],
	Sin[e2]*Sin[e1-e3]}
	]

Euler1232EP[e_] := Block[{c1,s1,c2,s2,c3,s3},
	c1 = Cos[e[[1]]/2];
	s1 = Sin[e[[1]]/2];
	c2 = Cos[e[[2]]/2];
	s2 = Sin[e[[2]]/2];
	c3 = Cos[e[[3]]/2];
	s3 = Sin[e[[3]]/2];
	
	{c1*c2*c3-s1*s2*s3,
	s1*c2*c3+c1*s2*s3,
	c1*s2*c3-s1*c2*s3,
	c1*c2*s3+s1*s2*c3}
	]

Euler1312EP[e_] := Block[{e1,e2,e3},
	e1 = e[[1]]/2;
	e2 = e[[2]]/2;
	e3 = e[[3]]/2;
	
	{Cos[e2]*Cos[e1+e3],
	Cos[e2]*Sin[e1+e3],
	Sin[e2]*Sin[-e1+e3],
	Sin[e2]*Cos[-e1+e3]}
	]

Euler1322EP[e_] := Block[{c1,s1,c2,s2,c3,s3},
	c1 = Cos[e[[1]]/2];
	s1 = Sin[e[[1]]/2];
	c2 = Cos[e[[2]]/2];
	s2 = Sin[e[[2]]/2];
	c3 = Cos[e[[3]]/2];
	s3 = Sin[e[[3]]/2];
	
	{c1*c2*c3+s1*s2*s3,
	s1*c2*c3-c1*s2*s3,
	c1*c2*s3-s1*s2*c3,
	c1*s2*c3+s1*c2*s3}
	]

Euler2122EP[e_] := Block[{e1,e2,e3},
	e1 = e[[1]]/2;
	e2 = e[[2]]/2;
	e3 = e[[3]]/2;
	
	{Cos[e2]*Cos[e1+e3],
	Sin[e2]*Cos[-e1+e3],
	Cos[e2]*Sin[e1+e3],
	Sin[e2]*Sin[-e1+e3]}
	]

Euler2132EP[e_] := Block[{c1,s1,c2,s2,c3,s3},
	c1 = Cos[e[[1]]/2];
	s1 = Sin[e[[1]]/2];
	c2 = Cos[e[[2]]/2];
	s2 = Sin[e[[2]]/2];
	c3 = Cos[e[[3]]/2];
	s3 = Sin[e[[3]]/2];
	
	{c1*c2*c3+s1*s2*s3,
	c1*s2*c3+s1*c2*s3,
	s1*c2*c3-c1*s2*s3,
	c1*c2*s3-s1*s2*c3}
	]

Euler2312EP[e_] := Block[{c1,s1,c2,s2,c3,s3},
	c1 = Cos[e[[1]]/2];
	s1 = Sin[e[[1]]/2];
	c2 = Cos[e[[2]]/2];
	s2 = Sin[e[[2]]/2];
	c3 = Cos[e[[3]]/2];
	s3 = Sin[e[[3]]/2];
	
	{c1*c2*c3-s1*s2*s3,
	c1*c2*s3+s1*s2*c3,
	s1*c2*c3+c1*s2*s3,
	c1*s2*c3-s1*c2*s3}
	]

Euler2322EP[e_] := Block[{e1,e2,e3},
	e1 = e[[1]]/2;
	e2 = e[[2]]/2;
	e3 = e[[3]]/2;
	
	{Cos[e2]*Cos[e1+e3],
	Sin[e2]*Sin[e1-e3],
	Cos[e2]*Sin[e1+e3],
	Sin[e2]*Cos[e1-e3]}
	]

Euler3122EP[e_] := Block[{c1,s1,c2,s2,c3,s3},
	c1 = Cos[e[[1]]/2];
	s1 = Sin[e[[1]]/2];
	c2 = Cos[e[[2]]/2];
	s2 = Sin[e[[2]]/2];
	c3 = Cos[e[[3]]/2];
	s3 = Sin[e[[3]]/2];
	
	{c1*c2*c3-s1*s2*s3,
	c1*s2*c3-s1*c2*s3,
	c1*c2*s3+s1*s2*c3,
	s1*c2*c3+c1*s2*s3}
	]

Euler3132EP[e_] := Block[{e1,e2,e3},
	e1 = e[[1]]/2;
	e2 = e[[2]]/2;
	e3 = e[[3]]/2;
	
	{Cos[e2]*Cos[e1+e3],
	Sin[e2]*Cos[e1-e3],
	Sin[e2]*Sin[e1-e3],
	Cos[e2]*Sin[e1+e3]}
	]

Euler3212EP[e_] := Block[{c1,s1,c2,s2,c3,s3},
	c1 = Cos[e[[1]]/2];
	s1 = Sin[e[[1]]/2];
	c2 = Cos[e[[2]]/2];
	s2 = Sin[e[[2]]/2];
	c3 = Cos[e[[3]]/2];
	s3 = Sin[e[[3]]/2];
	
	{c1*c2*c3+s1*s2*s3,
	c1*c2*s3-s1*s2*c3,
	c1*s2*c3+s1*c2*s3,
	s1*c2*c3-c1*s2*s3}
	]

Euler3232EP[e_] := Block[{e1,e2,e3},
	e1 = e[[1]]/2;
	e2 = e[[2]]/2;
	e3 = e[[3]]/2;
	
	{Cos[e2]*Cos[e1+e3],
	Sin[e2]*Sin[-e1+e3],
	Sin[e2]*Cos[-e1+e3],
	Cos[e2]*Sin[e1+e3]}
	]

Euler1212Gibbs[e_] := EP2Gibbs[Euler1212EP[e]]

Euler1232Gibbs[e_] := EP2Gibbs[Euler1232EP[e]]

Euler1312Gibbs[e_] := EP2Gibbs[Euler1312EP[e]]

Euler1322Gibbs[e_] := EP2Gibbs[Euler1322EP[e]]

Euler2122Gibbs[e_] := EP2Gibbs[Euler2122EP[e]]

Euler2132Gibbs[e_] := EP2Gibbs[Euler2132EP[e]]

Euler2312Gibbs[e_] := EP2Gibbs[Euler2312EP[e]]

Euler2322Gibbs[e_] := EP2Gibbs[Euler2322EP[e]]

Euler3122Gibbs[e_] := EP2Gibbs[Euler3122EP[e]]

Euler3132Gibbs[e_] := EP2Gibbs[Euler3132EP[e]]

Euler3232Gibbs[e_] := EP2Gibbs[Euler3232EP[e]]

Euler3212Gibbs[e_] := EP2Gibbs[Euler3212EP[e]]

Euler1212MRP[e_] := EP2MRP[Euler1212EP[e]]

Euler1232MRP[e_] := EP2MRP[Euler1232EP[e]]

Euler1312MRP[e_] := EP2MRP[Euler1312EP[e]]

Euler1322MRP[e_] := EP2MRP[Euler1322EP[e]]

Euler2122MRP[e_] := EP2MRP[Euler2122EP[e]]

Euler2132MRP[e_] := EP2MRP[Euler2132EP[e]]

Euler2312MRP[e_] := EP2MRP[Euler2312EP[e]]

Euler2322MRP[e_] := EP2MRP[Euler2322EP[e]]

Euler3122MRP[e_] := EP2MRP[Euler3122EP[e]]

Euler3132MRP[e_] := EP2MRP[Euler3132EP[e]]

Euler3232MRP[e_] := EP2MRP[Euler3232EP[e]]

Euler3212MRP[e_] := EP2MRP[Euler3212EP[e]]

Euler1212PRV[e_] := EP2PRV[Euler1212EP[e]]

Euler1232PRV[e_] := EP2PRV[Euler1232EP[e]]

Euler1312PRV[e_] := EP2PRV[Euler1312EP[e]]

Euler1322PRV[e_] := EP2PRV[Euler1322EP[e]]

Euler2122PRV[e_] := EP2PRV[Euler2122EP[e]]

Euler2132PRV[e_] := EP2PRV[Euler2132EP[e]]

Euler2312PRV[e_] := EP2PRV[Euler2312EP[e]]

Euler2322PRV[e_] := EP2PRV[Euler2322EP[e]]

Euler3122PRV[e_] := EP2PRV[Euler3122EP[e]]

Euler3132PRV[e_] := EP2PRV[Euler3132EP[e]]

Euler3232PRV[e_] := EP2PRV[Euler3232EP[e]]

Euler3212PRV[e_] := EP2PRV[Euler3212EP[e]]

Gibbs2C[q_] := Block[{q1,q2,q3,d1},
	q1 = q[[1]];
	q2 = q[[2]];
	q3 = q[[3]];
	d1 = Dot[q,q];
	
	{{1+2*q1*q1-d1,
	2*(q1*q2+q3),
	2*(q1*q3-q2)},{
	2*(q2*q1-q3),
	1+2*q2*q2-d1,
	2*(q2*q3+q1)},{
	2*(q3*q1+q2),
	2*(q3*q2-q1),
	1+2*q3*q3-d1}}/(1+d1)
	]

Gibbs2EP[q1_] := Block[{d},
	d = 1/Sqrt[1+Dot[q1,q1]];
	{d,
	q1[[1]]*d,q1[[2]]*d,q1[[3]]*d}
	]

Gibbs2Euler121[q_]:=EP2Euler121[Gibbs2EP[q]];

Gibbs2Euler123[q_]:=EP2Euler123[Gibbs2EP[q]];

Gibbs2Euler131[q_]:=EP2Euler131[Gibbs2EP[q]];

Gibbs2Euler132[q_]:=EP2Euler132[Gibbs2EP[q]];

Gibbs2Euler232[q_]:=EP2Euler232[Gibbs2EP[q]];

Gibbs2Euler231[q_]:=EP2Euler231[Gibbs2EP[q]];

Gibbs2Euler212[q_]:=EP2Euler212[Gibbs2EP[q]];

Gibbs2Euler213[q_]:=EP2Euler213[Gibbs2EP[q]];

Gibbs2Euler313[q_]:=EP2Euler313[Gibbs2EP[q]];

Gibbs2Euler312[q_]:=EP2Euler312[Gibbs2EP[q]];

Gibbs2Euler323[q_]:=EP2Euler323[Gibbs2EP[q]];

Gibbs2Euler321[q_]:=EP2Euler321[Gibbs2EP[q]];

Gibbs2MRP[q_]:=q/(1+Sqrt[1+Dot[q,q]]);

Gibbs2PRV[q1_]:= Block[{tp,p},
	tp = Sqrt[Dot[q1,q1]];
	p = 2*ArcTan[tp];
	q1/tp*p
	];

MRP2C[q_]:= Block[{q1,q2,q3,d1,S,e},
	q1 = q[[1]];
	q2 = q[[2]];
	q3 = q[[3]];
	
	d1 = Dot[q,q];
	S = 1-d1;
	d = (1+d1)*(1+d1);
	
	{{4*(2*q1*q1-d1)+S*S,
	8*q1*q2+4*q3*S,
	8*q1*q3-4*q2*S},{
	8*q2*q1-4*q3*S,
	4*(2*q2*q2-d1)+S*S,
	8*q2*q3+4*q1*S},{
	8*q3*q1+4*q2*S,
	8*q3*q2-4*q1*S,
	4*(2*q3*q3-d1)+S*S}}/d
	]

MRP2EP[q1_]:=Block[{ps},
	ps = 1+Dot[q1,q1];
	
	{(1-Dot[q1,q1])/ps,
	2*q1[[1]]/ps,
	2*q1[[2]]/ps,
	2*q1[[3]]/ps}
	]

MRP2Euler121[q_]:=EP2Euler121[MRP2EP[q]];

MRP2Euler123[q_]:=EP2Euler123[MRP2EP[q]];

MRP2Euler131[q_]:=EP2Euler131[MRP2EP[q]];

MRP2Euler132[q_]:=EP2Euler132[MRP2EP[q]];

MRP2Euler232[q_]:=EP2Euler232[MRP2EP[q]];

MRP2Euler231[q_]:=EP2Euler231[MRP2EP[q]];

MRP2Euler212[q_]:=EP2Euler212[MRP2EP[q]];

MRP2Euler213[q_]:=EP2Euler213[MRP2EP[q]];

MRP2Euler313[q_]:=EP2Euler313[MRP2EP[q]];

MRP2Euler312[q_]:=EP2Euler312[MRP2EP[q]];

MRP2Euler323[q_]:=EP2Euler323[MRP2EP[q]];

MRP2Euler321[q_]:=EP2Euler321[MRP2EP[q]];

MRP2Gibbs[q1_]:=2*q1/(1-Dot[q1,q1]);

MRP2PRV[q1_]:= Block[{tp,p},
	tp = Sqrt[Dot[q1,q1]];
	p = 4*ArcTan[tp];
	q1/tp*p
	];

MRPswitch[q_,s2_]:= Block[{q2},
	q2 = Dot[q,q];
	If[q2>s2*s2,-q/q2,q]
	]

Picheck[x_]:=Block[{q},
	q = x;
	If[x>Pi,q=x-2*Pi];
	If[x<-Pi,q=x+2*Pi];
	q
	]

PRV2C[q_]:=Block[{q0,q1,q2,q3,cp,sp,d1},
	q0 = Sqrt[Dot[q,q]];
	q1 = q[[1]]/q0;
	q2 = q[[2]]/q0;
	q3 = q[[3]]/q0;
	
	cp = Cos[q0];
	sp = Sin[q0];
	d1 = 1-cp;
	 
	{{q1*q1*d1+cp,
	q1*q2*d1+q3*sp,
	q1*q3*d1-q2*sp},{
	q2*q1*d1-q3*sp,
	q2*q2*d1+cp,
	q2*q3*d1+q1*sp},{
	q3*q1*d1+q2*sp,
	q3*q2*d1-q1*sp,
	q3*q3*d1+cp}}
	]
	
PRV2elem[r_]:=Block[{q1},
	q1 = Sqrt[Dot[r,r]];
	{q1,
	r[[1]]/q1,r[[2]]/q1,r[[3]]/q1}
	]

PRV2EP[q_]:=Block[{q1,sp},
	q1 = PRV2elem[q];
	sp = Sin[q1[[1]]/2];
	{Cos[q1[[1]]/2],
	q1[[2]]*sp,
	q1[[3]]*sp,
	q1[[4]]*sp}
	]

PRV2Euler121[q_]:= EP2Euler121[PRV2EP[q]]

PRV2Euler123[q_]:= EP2Euler123[PRV2EP[q]]

PRV2Euler131[q_]:= EP2Euler131[PRV2EP[q]]

PRV2Euler132[q_]:= EP2Euler132[PRV2EP[q]]

PRV2Euler232[q_]:= EP2Euler232[PRV2EP[q]]

PRV2Euler231[q_]:= EP2Euler231[PRV2EP[q]]

PRV2Euler212[q_]:= EP2Euler212[PRV2EP[q]]

PRV2Euler213[q_]:= EP2Euler213[PRV2EP[q]]

PRV2Euler313[q_]:= EP2Euler313[PRV2EP[q]]

PRV2Euler312[q_]:= EP2Euler312[PRV2EP[q]]

PRV2Euler323[q_]:= EP2Euler323[PRV2EP[q]]

PRV2Euler321[q_]:= EP2Euler321[PRV2EP[q]]

PRV2Gibbs[q_] := Block[{q1,tp},
	q1 = PRV2elem[q];
	tp = Tan[q1[[1]]/2];
	{q1[[2]]*tp,q1[[3]]*tp,q1[[4]]*tp}
	]
 
PRV2MRP[q_] := Block[{q1,tp},
	q1 = PRV2elem[q];
	tp = Tan[q1[[1]]/4];
	{q1[[2]]*tp,q1[[3]]*tp,q1[[4]]*tp}
	]

subEP[b1_,b2_] := {
	b2[[1]]*b1[[1]]+
	b2[[2]]*b1[[2]]+b2[[3]]*b1[[3]] + b2[[4]]*b1[[4]],
	-b2[[2]]*b1[[1]] + b2[[1]]*b1[[2]] +
	b2[[4]]*b1[[3]] - b2[[3]]*b1[[4]],
	-b2[[3]]*b1[[1]] - b2[[4]]*b1[[2]] +
	b2[[1]]*b1[[3]]+b2[[2]]*b1[[4]],
	-b2[[4]]*b1[[1]] + b2[[3]]*b1[[2]] -
	b2[[2]]*b1[[3]] + b2[[1]]*b1[[4]]}
	
subEuler121[e_,e1_] := Block[{cp,cp1,sp,sp1,
	dum,cp2,dum2},
	cp = Cos[e[[2]]];
	cp1 = Cos[e1[[2]]];
	sp = Sin[e[[2]]];
	sp1 = Sin[e1[[2]]];
	dum = e[[1]]-e1[[1]];
	dum2 = ArcCos[cp1*cp+sp1*sp*Cos[dum]];
	cp2 = Cos[dum2];
	{Picheck[-e1[[3]]+ArcTan[cp2*cp1-cp,
	sp1*sp*Sin[dum]]],
	dum2,
	Picheck[e[[3]]-ArcTan[cp1-cp*cp2,sp1*sp*Sin[dum]]]}
	]

subEuler131[e_,e1_] := Block[{cp,cp1,sp,sp1,
	dum,cp2,dum2},
	cp = Cos[e[[2]]];
	cp1 = Cos[e1[[2]]];
	sp = Sin[e[[2]]];
	sp1 = Sin[e1[[2]]];
	dum = e[[1]]-e1[[1]];
	dum2 = ArcCos[cp1*cp+sp1*sp*Cos[dum]];
	cp2 = Cos[dum2];
	{Picheck[-e1[[3]]+ArcTan[cp2*cp1-cp,
	sp1*sp*Sin[dum]]],
	dum2,
	Picheck[e[[3]]-ArcTan[cp1-cp*cp2,sp1*sp*Sin[dum]]]}
	]

subEuler212[e_,e1_] := Block[{cp,cp1,sp,sp1,
	dum,cp2,dum2},
	cp = Cos[e[[2]]];
	cp1 = Cos[e1[[2]]];
	sp = Sin[e[[2]]];
	sp1 = Sin[e1[[2]]];
	dum = e[[1]]-e1[[1]];
	dum2 = ArcCos[cp1*cp+sp1*sp*Cos[dum]];
	cp2 = Cos[dum2];
	{Picheck[-e1[[3]]+ArcTan[cp2*cp1-cp,
	sp1*sp*Sin[dum]]],
	dum2,
	Picheck[e[[3]]-ArcTan[cp1-cp*cp2,sp1*sp*Sin[dum]]]}
	]

subEuler232[e_,e1_] := Block[{cp,cp1,sp,sp1,
	dum,cp2,dum2},
	cp = Cos[e[[2]]];
	cp1 = Cos[e1[[2]]];
	sp = Sin[e[[2]]];
	sp1 = Sin[e1[[2]]];
	dum = e[[1]]-e1[[1]];
	dum2 = ArcCos[cp1*cp+sp1*sp*Cos[dum]];
	cp2 = Cos[dum2];
	{Picheck[-e1[[3]]+ArcTan[cp2*cp1-cp,
	sp1*sp*Sin[dum]]],
	dum2,
	Picheck[e[[3]]-ArcTan[cp1-cp*cp2,sp1*sp*Sin[dum]]]}
	]

subEuler313[e_,e1_] := Block[{cp,cp1,sp,sp1,
	dum,cp2,dum2},
	cp = Cos[e[[2]]];
	cp1 = Cos[e1[[2]]];
	sp = Sin[e[[2]]];
	sp1 = Sin[e1[[2]]];
	dum = e[[1]]-e1[[1]];
	dum2 = ArcCos[cp1*cp+sp1*sp*Cos[dum]];
	cp2 = Cos[dum2];
	{Picheck[-e1[[3]]+ArcTan[cp2*cp1-cp,
	sp1*sp*Sin[dum]]],
	dum2,
	Picheck[e[[3]]-ArcTan[cp1-cp*cp2,sp1*sp*Sin[dum]]]}
	]

subEuler323[e_,e1_] := Block[{cp,cp1,sp,sp1,
	dum,cp2,dum2},
	cp = Cos[e[[2]]];
	cp1 = Cos[e1[[2]]];
	sp = Sin[e[[2]]];
	sp1 = Sin[e1[[2]]];
	dum = e[[1]]-e1[[1]];
	dum2 = ArcCos[cp1*cp+sp1*sp*Cos[dum]];
	cp2 = Cos[dum2];
	{Picheck[-e1[[3]]+ArcTan[cp2*cp1-cp,
	sp1*sp*Sin[dum]]],
	dum2,
	Picheck[e[[3]]-ArcTan[cp1-cp*cp2,sp1*sp*Sin[dum]]]}
	]

subEuler123[e_,e1_] := Block[{C0,C1,C2},
	C0 = Euler1232C[e];
	C1 = Euler1232C[e1];
	C2 = C0.Transpose[C1];
	C2Euler123[C2]
	]

subEuler132[e_,e1_] := Block[{C0,C1,C2},
	C0 = Euler1322C[e];
	C1 = Euler1322C[e1];
	C2 = C0.Transpose[C1];
	C2Euler132[C2]
	]

subEuler213[e_,e1_] := Block[{C0,C1,C2},
	C0 = Euler2132C[e];
	C1 = Euler2132C[e1];
	C2 = C0.Transpose[C1];
	C2Euler213[C2]
	]

subEuler231[e_,e1_] := Block[{C0,C1,C2},
	C0 = Euler2312C[e];
	C1 = Euler2312C[e1];
	C2 = C0.Transpose[C1];
	C2Euler231[C2]
	]

subEuler312[e_,e1_] := Block[{C0,C1,C2},
	C0 = Euler3122C[e];
	C1 = Euler3122C[e1];
	C2 = C0.Transpose[C1];
	C2Euler312[C2]
	]

subEuler321[e_,e1_] := Block[{C0,C1,C2},
	C0 = Euler3212C[e];
	C1 = Euler3212C[e1];
	C2 = C0.Transpose[C1];
	C2Euler321[C2]
	]

subGibbs[q1_,q2_]:= (q1-q2+Cross[q1,q2])/(1+Dot[q1,q2])

subMRP[q1_,q2_] := ((1-Dot[q2,q2])*q1 - 
	(1-Dot[q1,q1])*q2 + 2*Cross[q1,q2])/
	(1+Dot[q1,q1]*Dot[q2,q2]+2*Dot[q1,q2])

subPRV[q11_,q22_]:= Block[{q1,q2,cp1,cp2,sp1,sp2,e1,e2,
	p,sp,e},
	
	q1 = PRV2elem[q11];
	q2 = PRV2elem[q22];
	cp1 = Cos[q1[[1]]/2];
	cp2 = Cos[q2[[1]]/2];
	sp1 = Sin[q1[[1]]/2];
	sp2 = Sin[q2[[1]]/2];
	e1 = {q1[[2]],q1[[3]],q1[[4]]};
	e2 = {q2[[2]],q2[[3]],q2[[4]]};
	
	p = 2*ArcCos[cp1*cp2+sp1*sp2*Dot[e1,e2]];
	sp = Sin[p/2];
	e = (-cp1*sp2*e2+cp2*sp1*e1+sp1*sp2
		*Cross[e1,e2])/sp;
	p*e
	]

tilde[q_] := {{0, -q[[3]], q[[2]]},{q[[3]],0,-q[[1]]},{-q[[2]],q[[1]],0}}


End[]

EndPackage[]
