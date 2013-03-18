/*
 *  orbitalMotion.c
 *  OrbitalMotion
 *
 *  Created by Hanspeter Schaub on 6/19/05.
 *  Copyright (c) 2005 Hanspeter Schaub. All rights reserved.
 *
 */

#include "orbitalMotion.h"



/*

%    f = E2f(Ecc, e)
%   
%   Maps eccentric anomaly angles into true anomaly angles.  
%   This function requires the orbit to be either circular or 
%   non-rectilinar elliptic orbit.
%
%   Input is
%       Ecc - eccentric anomaly    (rad)
%       e   - eccentricity          0 <= e < 1
%
%   Output is
%       f   - true anomaly         (rad)

*/
double E2f(double Ecc, double e) 
{
    double f;
    
    if ((e >= 0) && (e < 1))
        f = 2*atan2(sqrt(1+e)*sin(Ecc/2), sqrt(1-e)*cos(Ecc/2));    
    else {
        f = NAN;
        printf("ERROR: E2f() received e = %g \n",e);
        printf("The value of e should be 0 <= e < 1. \n");
    }
    
    return f;
}









/*
%   M = E2M(Ecc, e)
%   
%   Maps the eccentric anomaly angle into the corresponding
%   mean elliptic anomaly angle.  Both 2D and 1D elliptic
%   orbit are allowed.  
%
%   Input is
%       Ecc - eccentric anomaly    (rad)
%       e   - eccentricity          0 <= e <= 1
%
%   Output is
%       M   - mean elliptic anomaly (rad)

*/
double E2M(double Ecc, double e) 
{
    double M;
    
    if ((e >= 0) && (e < 1))
        M = Ecc - e*sin(Ecc);    
    else {
        M = NAN;
        printf("ERROR: E2M() received e = %g \n",e);
        printf("The value of e should be 0 <= e < 1. \n");
    }
    
    return M;
}









/*
%   Ecc = f2E(f, e)
%   
%   Maps true anomaly angles into eccentric anomaly angles.  
%   This function requires the orbit to be either circular or 
%   non-rectilinar elliptic orbit.
%
%   Input is
%       f   - true anomaly angle   (rad)
%       e   - eccentricity          0 <= e < 1
%
%   Output is
%       Ecc - eccentric anomaly     (rad)

*/
double f2E(double f, double e) 
{
    double Ecc;
    
    if ((e >= 0) && (e < 1))
        Ecc = 2*atan2(sqrt(1-e)*sin(f/2), sqrt(1+e)*cos(f/2));    
    else {
        Ecc = NAN;
        printf("ERROR: f2E() received e = %g \n",e);
        printf("The value of e should be 0 <= e < 1. \n");
    }
    
    return Ecc;
}









/*
%   H = f2H(f, e)
%   
%   Maps true anomaly angles into hyperbolic anomaly angles.  
%   This function requires the orbit to be hyperbolic
%
%   Input is
%       f   - true anomaly angle   (rad)
%       e   - eccentricity          e > 1
%
%   Output is
%       H   - hyperbolic anomaly   (rad)

*/
double f2H(double f, double e) 
{
    double H;
    
    if (e > 1)
        H = 2*arc_tanh(sqrt((e-1)/(e+1))*tan(f/2));    
    else {
        H = NAN;
        printf("ERROR: f2H() received e = %g \n",e);
        printf("The value of e should be 1 < e. \n");
    }
    
    return H;
}









/*
%   f = H2f(H, e)
%   
%   Maps hyperbolic anomaly angles into true anomaly angles.  
%   This function requires the orbit to be hyperbolic
%
%   Input is
%       H   - hyperbolic anomaly   (rad)
%       e   - eccentricity          e > 1
%
%   Output is
%       f   - true anomaly         (rad)

*/
double H2f(double H, double e) 
{
    double f;
    
    if (e > 1)
        f = 2*atan(sqrt((e+1)/(e-1))*tanh(H/2));    
    else {
        f = NAN;
        printf("ERROR: H2f() received e = %g \n",e);
        printf("The value of e should be 1 < e \n");
    }
    
    return f;
}









/*
%   N = H2N(H, e)
%   
%   Maps the hyperbolic anomaly angle H into the corresponding
%   mean hyperbolic anomaly angle N.  
%
%   Input is
%       H   - hyperbolic anomaly        (rad)
%       e   - eccentricity              e > 1
%
%   Output is
%       N   - mean hyperbolic anomaly   (rad)
*/
double H2N(double H, double e) 
{
    double N;
    
    if (e > 1)
        N = e*sinh(H) - H;    
    else {
        N = NAN;
        printf("ERROR: H2N() received e = %g \n",e);
        printf("The value of e should be 1 < e. \n");
    }
    
    return N;
}









/*
%   Ecc = M2E(M, e)
%   
%   Maps the mean elliptic anomaly angle into the corresponding
%   eccentric anomaly angle.  Both 2D and 1D elliptic
%   orbit are allowed.  
%
%   Input is
%       M   - mean elliptic anomaly     (rad)
%       e   - eccentricity              0 <= e <= 1
%
%   Output is
%       Ecc - eccentric anomaly         (rad)
*/
double M2E(double M, double e) 
{
    double small, dE,E1;
	int max = 200,
		count = 0;
	small = 1e-13;
    dE = 10*small;
	E1 = M;
	
    if ((e >= 0) && (e < 1))
        while(fabs(dE)>small) {
            dE = (E1-e*sin(E1)-M)/(1-e*cos(E1));
            E1 -= dE;
            if (++count > max) {
                printf("iteration error in M2E(%f,%f)\n",M,e);
                dE = 0.;
            }
        }    
    else {
        E1 = NAN;
        printf("ERROR: M2E() received e = %g \n",e);
        printf("The value of e should be 0 <= e <= 1. \n");
    }
    
    return E1;
}









/*
%   H = N2H(N, e)
%   
%   Maps the mean hyperbolic anomaly angle N into the corresponding
%   hyperbolic anomaly angle H.  
%
%   Input is
%       N   - mean hyperbolic anomaly   (rad)
%       e   - eccentricity              e > 1
%
%   Output is
%       H   - hyperbolic anomaly        (rad)
*/
double N2H(double N, double e) 
{
    double small, dH,H1;
	int max = 200,
		count = 0;
	small = 1e-13;
    dH = 10*small;
	H1 = N;
	
    if (e > 1)
        while(fabs(dH)>small) {
            dH = (e*sinh(H1) - H1 - N)/(e*cosh(H1)-1);
            H1 -= dH;
            if (++count > max) {
                printf("iteration error in N2H(%f,%f)\n",N,e);
                dH = 0.;
            }
        }    
    else {
        H1 = NAN;
        printf("ERROR: N2H() received e = %g \n",e);
        printf("The value of e should be 0 <= e <= 1. \n");
    }
    
    return H1;
}








/*
%   elem2rv(mu, *element, *rVec, *vVec)
%   
%   Translates the orbit elements
%            a   - semi-major axis           (km)
%            e   - eccentricity
%            i   - inclination               (rad)
%            AN  - ascending node            (rad)
%            AP  - argument of periapses     (rad)
%            f   - true anomaly angle        (rad)
%   to the inertial Cartesian position and velocity vectors.
%   The attracting body is specified through the supplied
%   gravitational constant mu (units of km^3/s^2).
%
%   The code can handle the following cases:
%       circular:       e = 0           a > 0
%       elliptical-2D:  0 < e < 1       a > 0
%       elliptical-1D:  e = 1           a > 0        f = Ecc. Anom. here
%       parabolic:      e = 1           rp = -a
%       hyperbolic:     e > 1           a < 0
%
%   Note: to handle the parabolic case and distinguish it form the
%   rectilinear elliptical case, instead of passing along the
%   semi-major axis a in the "a" input slot, the negative radius
%   at periapses is supplied.  Having "a" be negative and e = 1
%   is a then a unique identified for the code for the parabolic 
%   case.
*/
void elem2rv(double mu, classicElements *elements, double *rVec, double *vVec)
{
    double e, a, Ecc, f, r, v, i, rp, p, AP, AN, theta, h;
    double ir[3+1];
    
    /* map classical elements structure into local variables */
    a  = elements->a;
    e  = elements->e;
    i  = elements->i;
    AN = elements->Omega;
    AP = elements->omega;
    f  = elements->anom;
    
    if ((e == 1) && (a > 0)) {    /* rectilinear elliptic orbit case */
        Ecc = f;                /* f is treated as ecc. anomaly */
        r = a*(1-e*cos(Ecc));   /* orbit radius  */
        v = sqrt(2*mu/r - mu/a);
        ir[1] = cos(AN)*cos(AP) - sin(AN)*sin(AP)*cos(i);
        ir[2] = sin(AN)*cos(AP) + cos(AN)*sin(AP)*cos(i);
        ir[3] = sin(AP)*sin(i);
        mult(r,ir, rVec);
        if (sin(Ecc)>0)
            mult(-v,ir,vVec);
        else
            mult(v,ir,vVec);
        
        
    } else {

        if ((e == 1) && (a < 0)) {   /* parabolic case */
            rp = -a;                /* radius at periapses  */
            p = 2*rp;               /* semi-latus rectum */
        } else {                    /* elliptic and hyperbolic cases */
            p = a*(1-e*e);    /* semi-latus rectum */
        }
        
        r       = p/(1+e*cos(f));   /* orbit radius */
        theta   = AP + f;           /* true latitude angle */
        h       = sqrt(mu*p);       /* orbit ang. momentum mag. */
                
        rVec[1] = r*(cos(AN)*cos(theta) - sin(AN)*sin(theta)*cos(i));
        rVec[2] = r*(sin(AN)*cos(theta) + cos(AN)*sin(theta)*cos(i));
        rVec[3] = r*(sin(theta)*sin(i));
        
        vVec[1] =-mu/h*(cos(AN)*(sin(theta)+e*sin(AP)) + sin(AN)*(cos(theta)+e*cos(AP))*cos(i));
        vVec[2] =-mu/h*(sin(AN)*(sin(theta)+e*sin(AP)) - cos(AN)*(cos(theta)+e*cos(AP))*cos(i));
        vVec[3] =-mu/h*(-(cos(theta)+e*cos(AP))*sin(i));

    }
    
    return;
}








/*
%   rv2elem(mu, *rVec, *vVec, *elements)
%   
%   Translates the orbit elements inertial Cartesian position 
%   vector rVec and velocity vector vVec into the corresponding 
%   classical orbit elements where
%            a   - semi-major axis           (km)
%            e   - eccentricity
%            i   - inclination               (rad)
%            AN  - ascending node            (rad)
%            AP  - argument of periapses     (rad)
%            f   - true anomaly angle        (rad)  
%                  if the orbit is rectilinear, then this will be the
%                  eccentric or hyperbolic anomaly
%   The attracting body is specified through the supplied
%   gravitational constant mu (units of km^3/s^2).
%
%   The code can handle the following cases:
%       circular:       e = 0           a > 0
%       elliptical-2D:  0 < e < 1       a > 0
%       elliptical-1D:  e = 1           a > 0
%       parabolic:      e = 1           a = -rp
%       hyperbolic:     e > 1           a < 0
%
%   For the parabolic case the semi-major axis is not defined.  
%   In this case -rp (radius at periapses) is returned instead
%   of a.  For the circular case, the AN and AP are ill-defined,
%   along with the associated ie and ip unit direction vectors
%   of the perifocal frame. In this circular orbit case, the 
%   unit vector ie is set equal to the normalized inertial 
%   position vector ir.  
*/
void rv2elem(double mu, double *rVec, double *vVec, classicElements *elements)
{
    double  r, h, eps, ai, p, rp, Ecc, H;
    double  ir[3+1], hVec[3+1], cVec[3+1], dum[3+1], ih[3+1], 
            ie[3+1], ip[3+1], dum2[3+1];
    
    /* define a small number */
    eps = 0.000000000001;

    /* compute orbit radius */
    r = norm(rVec);
    mult(1./r,rVec, ir);
    
    /* compute the angular momentum vector */
    cross(rVec,vVec, hVec);
    h    = norm(hVec);

    /* compute the eccentricity vector */
    cross(vVec, hVec, cVec);
    mult(- mu/r,rVec, dum);
    add(cVec, dum, cVec);
    elements->e    = norm(cVec)/mu;

    /* compute semi-major axis */
    ai = 2./r - dot(vVec,vVec)/mu;
    if (fabs(ai) > eps) {
        /* elliptic or hyperbolic case */
        elements->a = 1/ai;
    } else {
        /* parabolic case */
        p  = h*h/mu;
        rp = p/2;
        elements->a  = -rp;   /* a is not defined for parabola, so -rp is returned instead */
        elements->e  = 1;
    }
    
    if (h < eps) {  /* rectilinear motion case */
        equal(ir,ie);
        /* ip and ih are arbitrary */
        set(0,0,1,dum);
        set(0,1,0,dum2);
        cross(ie, dum,  ih);
        cross(ie, dum2, ip);
        if (norm(ih) > norm(ip)) {
            mult(1/norm(ih), ih, ih);
        } else {
            mult(1/norm(ip), ip, ih);
        }
        cross(ih, ie, ip);
        
    } else {
        /* compute perifocal frame unit direction vectors */
        mult(1./h, hVec, ih);
        if (fabs(elements->e) > eps) {
            /* non-circular case */
            mult(1./mu/elements->e, cVec, ie);
        } else {
            /* circular orbit case.  Here ie, ip are arbitrary, as long as they
               are perpenticular to the ih vector.  */
            equal(ir,ie);
        }
        cross(ih,ie,ip);
    }
    
    /* compute the 3-1-3 orbit plane orientation angles */
    elements->Omega = atan2(ih[1],-ih[2]);
    elements->i     = acos(ih[3]);
    elements->omega = atan2(ie[3],ip[3]);

    if (h < eps) {          /* rectilinear motion case */
        if (ai > 0) {       /* elliptic case */
            Ecc = acos(1-r*ai);
            if (dot(rVec,vVec) > 0) 
                Ecc = 2*M_PI - Ecc;
            elements->anom = Ecc;    /* for this mode the eccentric anomaly is returned */
        } else {            /* hyperbolic case */
            H = arc_cosh(r*ai+1);
            if (dot(rVec,vVec) < 0) 
                H = 2*M_PI - H;
            elements->anom = H;    /* for this mode the hyperbolic anomaly is returned */
        }
    } else {
        /* compute true anomaly */
        cross(ie,ir,dum);
        elements->anom = atan2(dot(dum,ih),dot(ie,ir));   
    }
    
    return;
}













/*
%   AtmosphericDensity (double alt)
%       
%   Purpose:   This program computes the atmospheric density based on altitude 
%              supplied by user.  This function uses a curve fit based on
%              atmospheric data from the Standard Atmoshere 1976 Data. This 
%              function is valid for altitudes ranging from 100km to 1000km. 
%
%              Note: This code can only be applied to spacecraft orbiting 
%              the Earth
%
%   Curve fit equation based on data obtained from: 
%   U.S. Standard Atmosphere, 1976, U.S. Government Printing Office, Washington, D.C., 1976.
%
%   Input is
%       alt -  This is the altitude supplied by the user in km
%
%   Output is
%       density  - This is the density at the given altitude in kg/m^3
*/
double AtmosphericDensity( double alt)
{
    double logdensity, density, val;
    
    /* Smooth exponential drop-off after 1000 km */
    if (alt > 1000.) {
        logdensity = (-7e-05)*alt-14.464;
        density = pow(10.,logdensity);
        return density;
    }
     
    /* Calculating the density based on a scaled 6th order polynomial fit 
       to the log of density */
    val = (alt-526.8000)/292.8563;
    logdensity = 0.34047*pow(val,6) - 0.5889*pow(val,5) -0.5269*pow(val,4) 
                + 1.0036*pow(val,3) + 0.60713*pow(val,2) -2.3024*val -12.575;
    
    /* Calculating density by raising 10 to the log of density */
    density = pow(10.,logdensity);
    
    return density;
}







/*
%   AtmosphericDrag (Cd,A,m,rvec,vvec, advec)
%       
%   Purpose:   This program computes the atmospheric drag acceleration
%              vector acting on a spacecraft.
%              Note the acceleration vector output is inertial, and is 
%              only valid for altitudes up to 1000 km.
%              Afterwards the drag force is zero. Only valid for Earth.
%
%   Input is
%       Cd -  This is the drag coefficient of the spacecraft
%       A  -  This is the cross-sectional area of the spacecraft in m^2
%       m  -  This is the mass of the spacecraft in kg
%       rvec - Inertial position vector of the spacecraft in km  [x;y;z]
%       vvec - Inertial velocity vector of the spacecraft in km/s [vx;vy;vz]
%
%   Output is
%       advec  - The inertial acceleration vector due to atmospheric 
%                drag in km/sec^2
%
*/
void   AtmosphericDrag(double Cd, double A, double m, double *rvec, 
                       double *vvec, double *advec)
{
    double r, v, alt, ad, density;
    
    /* find the altitude and velocity */
    r   = norm(rvec);
    v   = norm(vvec);
    alt = r - REQ_EARTH;

    /* %Checking if user supplied a orbital position is inside the earth */
    if (alt <= 0.) {
        printf("ERROR: AtmosphericDrag() received rvec = [%g %g %g] \n",rvec[1], rvec[2], rvec[3]);
        printf("The value of rvec should produce a positive altitude for the Earth.\n");
        set(NAN,NAN,NAN,advec);
        return;
    }
    
    /* get the Atmospheric density at the given altitude in kg/m^3 */
    density = AtmosphericDensity(alt);
    
    /* compute the magnitude of the drag acceleration */
    ad=((-0.5)*density*(Cd*A/m)*(pow(v*1000.,2)))/1000.;
    
    /* computing the vector for drag acceleration */
    mult(ad/v,vvec, advec);
    
    
    return;
}




/*
%   JPerturb(rvec,num,ajtot)
%       
%   Purpose:  Computes the J2-J6 zonal graviational perturbation
%             accelerations.  
%
%   Input is
%       rvec - Cartesian Position vector in kilometers [x;y;z].
%       num  - Corresponds to which J components to use,
%              must be an integer between 2 and 6. 
%              (note: Additive- 2 corresponds to J2 while 3 will
%              correspond to J2 + J3)
%   
%
%   Output is
%       ajtot  - The total acceleration vector due to the J 
%                perturbations in km/sec^2 [accelx;accely;accelz]
*/
void JPerturb(double *rvec, int num, double *ajtot)
{
    double mu, req, x, y, z, r, temp[4], temp2[4];
    
    /* Constants for Earth */
    mu  = MU_EARTH;
    req = REQ_EARTH;
        


    /* Calculate the J perturbations */
    x=rvec[1];
    y=rvec[2];
    z=rvec[3];
    r=norm(rvec);

    /* Error Checking */
    if ((num < 2) || (num > 6)) {
        printf("ERROR: JPerturb() received num = %d \n",num);
        printf("The value of num should be 2 <= num <= 6. \n");
        set(NAN,NAN,NAN,ajtot);
        return;
    }
    
    /* Calculating the total acceleration based on user input */
    if (num >= 2) {
        set((1-5*pow(z/r,2))*(x/r), (1-5*pow(z/r,2))*(y/r), (3-5*pow(z/r,2))*(z/r), ajtot);
        mult(-3./2.*J2*(mu/pow(r,2))*pow(req/r,2), ajtot, ajtot);
    } 
    if (num >= 3) {
        set(5*(7*pow(z/r,3)-3*(z/r))*(x/r), 5*(7*pow(z/r,3)-3*(z/r))*(y/r), -3*(10*pow(z/r,2)-(35./3.)*pow(z/r,4)-1), temp);
        mult(1./2.*J3*(mu/pow(r,2))*pow(req/r,3), temp, temp2);
        add(ajtot, temp2, ajtot); 
    }
    if (num >= 4) {
        set((3-42*pow(z/r,2)+63*pow(z/r,4))*(x/r), (3-42*pow(z/r,2)+63*pow(z/r,4))*(y/r), (15-70*pow(z/r,2)+63*pow(z/r,4))*(z/r), temp);
        mult(5./8.*J4*(mu/pow(r,2))*pow(req/r,4), temp, temp2);
        add(ajtot, temp2, ajtot);
    }
    if (num >= 5) {
        set(3*(35*(z/r)-210*pow(z/r,3)+231*pow(z/r,5))*(x/r), 3*(35*(z/r)-210*pow(z/r,3)+231*pow(z/r,5))*(y/r), -(15-315*pow(z/r,2)+945*pow(z/r,4)-693*pow(z/r,6)), temp);
        mult(1./8.*J5*(mu/pow(r,2))*pow(req/r,5), temp, temp2);
        add(ajtot, temp2, ajtot);
    }
    if (num >= 6) {
        set((35-945*pow(z/r,2)+3465*pow(z/r,4)-3003*pow(z/r,6))*(x/r), (35-945*pow(z/r,2)+3465*pow(z/r,4)-3003*pow(z/r,6))*(y/r), -(3003*pow(z/r,6)-4851*pow(z/r,4)+2205*pow(z/r,2)-245)*(z/r), temp);
        mult(-1./16.*J6*(mu/pow(r,2))*pow(req/r,6), temp, temp2);
        add(ajtot, temp2, ajtot);
    }


    return;
}



