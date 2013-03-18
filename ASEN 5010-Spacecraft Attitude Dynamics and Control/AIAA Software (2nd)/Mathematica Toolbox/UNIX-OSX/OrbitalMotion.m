(* :Title: OrbitalMotion *)

(* :Author: Hanspeter Schaub *)

(* :Summary:
This package provides various orbital mechanics subroutines using in astrodynamics calculations.

	
*)


(* :Package Version: 1.0 *)

(* :Copyright: Freeware *)

(* :History: 
	Version 1.0 by Hanspeter Schaub (Virginia Tech), June 2005.
*)

(* :Keywords: orbital motion, astrodynamics *)

(* :Mathematica Version: 5.0 *)

Needs["RigidBodyKinematics`"]
Needs["astroConstants`"]


BeginPackage["OrbitalMotion`"]

(* Usage messages *)

E2f::usage = "
    f = E2f[Ecc, e]\n
\n
Maps eccentric anomaly angles into true anomaly angles.  This function requires the orbit to be either circular or non-rectilinar elliptic orbit.\n
\n
   Input is\n
       Ecc - eccentric anomaly    (rad)\n
       e   - eccentricity          0 <= e < 1\n
\n
   Output is\n
       f   - true anomaly         (rad)\n"



E2M::usage = "
    M = E2M[Ecc, e]\n
   \n
Maps the eccentric anomaly angle into the corresponding mean elliptic anomaly angle.  Both 2D and 1D elliptic orbit are allowed.  \n
\n
    Input is\n
       Ecc - eccentric anomaly    (rad)\n
       e   - eccentricity          0 <= e <= 1\n
\n
   Output is\n
       M   - mean elliptic anomaly (rad)\n
"


f2E::usage = "
   Ecc = f2E[f, e]\n
  \n
   Maps true anomaly angles into eccentric anomaly angles.  This function requires the orbit to be either circular or non-rectilinar elliptic orbit.\n
\n
   Input is\n
       f   - true anomaly angle   (rad)\n
       e   - eccentricity          0 <= e < 1\n
\n
   Output is\n
       Ecc - eccentric anomaly     (rad)\n
"


f2H::usage = "
   H = f2H[f, e]\n
   \n
   Maps true anomaly angles into hyperbolic anomaly angles.  This function requires the orbit to be hyperbolic\n
\n
   Input is\n
       f   - true anomaly angle   (rad)\n
       e   - eccentricity          e > 1\n
\n
   Output is\n
       H   - hyperbolic anomaly   (rad)\n
"


H2f::usage = "
   f = H2f[H, e]\n
   \n
   Maps hyperbolic anomaly angles into true anomaly angles.  This function requires the orbit to be hyperbolic\n
\n
   Input is\n
       H   - hyperbolic anomaly   (rad)\n
       e   - eccentricity          e > 1\n
\n
   Output is\n
       f   - true anomaly         (rad)\n
"


H2N::usage = "
   N = H2N[H, e]\n
   \n
   Maps the hyperbolic anomaly angle H into the corresponding mean hyperbolic anomaly angle N.  \n
\n
   Input is\n
       H   - hyperbolic anomaly        (rad)\n
       e   - eccentricity              e > 1\n
\n
   Output is\n
       N   - mean hyperbolic anomaly   (rad)\n
"


M2E::usage = "
   Ecc = M2E[M, e]\n
   \n
   Maps the mean elliptic anomaly angle into the corresponding eccentric anomaly angle.  Both 2D and 1D elliptic orbit are allowed.  \n
\n
   Input is\n
       M   - mean elliptic anomaly     (rad)\n
       e   - eccentricity              0 <= e <= 1\n
\n
   Output is\n
       Ecc - eccentric anomaly         (rad)\n
"


N2H::usage = "
   H = N2H[N, e]\n
   \n
   Maps the mean hyperbolic anomaly angle N into the corresponding hyperbolic anomaly angle H.  \n
\n
   Input is\n
       N   - mean hyperbolic anomaly   (rad)\n
       e   - eccentricity              e > 1\n
\n
   Output is\n
       H   - hyperbolic anomaly        (rad)\n
"


elem2rv::usage = "
   {rVec, vVec} = elem2rv[mu, a, e, i, AN, AP, f]\n
   \n
   Translates the orbit elements\n
            a   - semi-major axis           (km)\n
            e   - eccentricity\n
            i   - inclination               (rad)\n
            AN  - ascending node            (rad)\n
            AP  - argument of periapses     (rad)\n
            f   - true anomaly angle        (rad)\n
   to the inertial Cartesian position and velocity vectors. The attracting body is specified through the supplied gravitational constant mu (units of km^3/s^2).\n
\n
   The code can handle the following cases:\n
       circular:       e = 0           a > 0\n
       elliptical-2D:  0 < e < 1       a > 0\n
       elliptical-1D:  e = 1           a > 0        f = Ecc. Anom. here\n
       parabolic:      e = 1           rp = -a\n
       hyperbolic:     e > 1           a < 0\n
   Note: to handle the parabolic case and distinguish it form the rectilinear elliptical case, instead of passing along the semi-major axis a in the 2nd input slot, the negative radius at periapses is supplied.  Having a be negative and e = 1 is a then a unique identified for the code for the parabolic case.\n
"

rv2elem::usage = "
   {a, e, i, AN, AP, f} = rv2elem[mu, rVec, vVec]\n
   \n
   Translates the orbit elements inertial Cartesian position vector rVec and velocity vector vVec into the corresponding classical orbit elements where\n
            a   - semi-major axis           (km)\n
            e   - eccentricity\n
            i   - inclination               (rad)\n
            AN  - ascending node            (rad)\n
            AP  - argument of periapses     (rad)\n
            f   - true anomaly angle        (rad)\n
   The attracting body is specified through the supplied gravitational constant mu (units of km^3/s^2).\n
\n
   The code can handle the following cases:\n
       circular:       e = 0           a > 0\n
       elliptical-2D:  0 < e < 1       a > 0\n
       elliptical-1D:  e = 1           a > 0\n
       parabolic:      e = 1           a = -rp\n
       hyperbolic:     e > 1           a < 0\n
\n
   For the parabolic case the semi-major axis is not defined.  In this case -rp (radius at periapses) is returned instead of a.  For the circular case, the AN and AP are ill-defined, along with the associated ie and ip unit direction vectors of the perifocal frame. In this circular orbit case, the unit vector ie is set equal to the normalized inertial position vector ir.  \n
"

AtmosphericDensity::usage = "
density = AtmosphericDensity[alt]\n
       \n
Purpose:   This program computes the atmospheric density based on altitude supplied by user.  This function uses a curve fit based on atmospheric data from the Standard Atmoshere 1976 Data. This function is valid for altitudes ranging from 100km to 1000km. \n
\n
Note: This code can only be applied to spacecraft orbiting the Earth\n
\n
Input is\n
   alt -  This is the altitude supplied by the user in km \n
\n
Output is\n
   density  - This is the density at the given altitude in kg/m^3\n
"




AtmosphericDrag::usage = "
{advec} = AtmosphericDrag [Cd,A,m,rvec,vvec]\n
  \n
Purpose:   This program computes the atmospheric drag acceleration vector acting on a spacecraft.  Note the acceleration vector output is inertial, and is only valid for altitudes up to 1000 km.  Afterwards the drag force is zero. Only valid for Earth.\n
\n
Input is\n
   Cd -  This is the drag coefficient of the spacecraft \n
   A  -  This is the cross-sectional area of the spacecraft in m^2 \n
   m  -  This is the mass of the spacecraft in kg \n
   rvec - Inertial position vector of the spacecraft in km  [x;y;z] \n
   vvec - Inertial velocity vector of the spacecraft in km/s [vx;vy;vz] \n
\n
Output is \n
   advec  - The inertial acceleration vector due to atmospheric \n
            drag in km/sec^2 \n
"


JPerturb::usage = "
{ajtot} = JPerturb[rvec,num] \n
 \n
Purpose:  Computes the J2-J6 zonal graviational perturbation accelerations. \n
 \n
Input is \n
    rvec - Cartesian Position vector in kilometers [x;y;z]. \n
    num  - Corresponds to which J components to use, must be an integer between 2 and 6. \n
(note: Additive- 2 corresponds to J2 while 3 will correspond to J2 + J3) \n
 \n
 \n
Output is \n
    ajtot  - The total acceleration vector due to the J perturbations in km/sec^2 {accelx;accely;accelz} \n
"










Begin["`Private`"]


E2f::"eccentricityCheck" = "Eccentricity `1` is out of bounds of [0,1)."
E2f[Ecc_, e_] := Block[{f},
    If [(e >= 0) && (e < 1), 
        f = 2*ArcTan[Sqrt[1-e]*Cos[Ecc/2],Sqrt[1+e]*Sin[Ecc/2]];
    ,
        f = Null;
        Message[E2f::"eccentricityCheck", e];
    ];
    Return[f];
]


E2M::"eccentricityCheck" = "Eccentricity `1` is out of bounds of [0,1)."
E2M[Ecc_, e_] := Block[{M},
    If [(e >= 0) && (e < 1), 
        M = Ecc - e*Sin[Ecc];
    ,
        M = Null;
        Message[E2M::"eccentricityCheck", e];
    ];
    Return[M];
]


f2E::"eccentricityCheck" = "Eccentricity `1` is out of bounds of [0,1)."
f2E[f_, e_] := Block[{Ecc},
    If [(e >= 0) && (e < 1), 
        Ecc = 2*ArcTan[Sqrt[1+e]*Cos[f/2],Sqrt[1-e]*Sin[f/2]];
    ,
        Ecc = Null;
        Message[f2E::"eccentricityCheck", e];
    ];
    Return[Ecc];
]


f2H::"eccentricityCheck" = "Eccentricity `1` is out of bounds of (1,infty)."
f2H[f_, e_] := Block[{H},
    If [e > 1, 
        H =  2*ArcTanh[Sqrt[(e-1)/(e+1)]*Tan[f/2]];
    ,
        H = Null;
        Message[f2H::"eccentricityCheck", e];
    ];
    Return[H];
]


H2f::"eccentricityCheck" = "Eccentricity `1` is out of bounds of (1,infty)."
H2f[H_, e_] := Block[{f},
    If [e > 1, 
        f =  2*ArcTan[Sqrt[(e+1)/(e-1)]*Tanh[H/2]];
    ,
        f = Null;
        Message[H2f::"eccentricityCheck", e];
    ];
    Return[f];
]


H2N::"eccentricityCheck" = "Eccentricity `1` is out of bounds of (1,infty)."
H2N[H_, e_] := Block[{NN},
    If [e > 1, 
        NN = e*Sinh[H] - H;
    ,
        NN = Null;
        Message[H2N::"eccentricityCheck", e];
    ];
    Return[NN];
]


M2E::"eccentricityCheck" = "Eccentricity `1` is out of bounds of [0,1]."
M2E[M_, e_] := Block[{Ecc},
    If [(e >= 0) && (e <= 1), 
        Ecc = (E2 /. FindRoot[M - (E2 - e Sin[E2]) == 0, {E2, M}]);
    ,
        Ecc = Null;
        Message[M2E::"eccentricityCheck", e];
    ];
    Return[Ecc];
]


N2H::"eccentricityCheck" = "Eccentricity `1` is out of bounds of (1,infty)."
N2H[NN_, e_] := Block[{H},
    If [e > 1, 
        H = (H2 /. FindRoot[NN - (e Sinh[H2] - H2) == 0, {H2, NN}]);
    ,
        H = Null;
        Message[N2H::"eccentricityCheck", e];
    ];
    Return[H];
]


elem2rv[mu_, a_, e_, i_, AN_, AP_, f_] := Block[{rp,p,theta,rVec,vVec},
    
    cAN = Cos[AN];      sAN = Sin[AN];
    ci  = Cos[i];       si  = Sin[i];
    cAP = Cos[AP];      sAP = Sin[AP];
    
    If[ (e == 1) && (a > 0),    (* rectilinear elliptic orbit case *)
        Ecc = f;                (* f is treated as ecc. anomaly *)
        r = a*(1-e*Cos[Ecc]);   (* orbit radius *)
        v = Sqrt[2*mu/r - mu/a];
        ir = {cAN*cAP - sAN*sAP*ci, sAN*cAP + cAN*sAP*ci, sAP*si};
        rVec = r*ir;
        If[Sin[Ecc]>0,
            vVec = -v*ir;
        ,
            vVec =  v*ir;
        ];
    ,
    
    If[(e == 1) && (a < 0), (* parabolic case *)
        rp = -a;            (* radius at periapses *)
        p = 2*rp;           (* semi-latus rectum *)
    ,                       (* elliptic and hyperbolic cases *)
        p = a*(1-e^2);      (* semi-latus rectum *)
    ];
    
    r = p/(1+e*Cos[f]);         (* orbit radius *)
    theta   = AP + f;           (* true latitude angle *)
    h       = Sqrt[mu*p];       (* orbit ang. momentum mag. *)
    
    cth = Cos[theta];   sth = Sin[theta];
    
    rVec    = r*{cAN*cth - sAN*sth*ci,
                 sAN*cth + cAN*sth*ci,
                 sth*si};
    
    vVec    =-mu/h*{cAN*(sth+e*sAP) + sAN*(cth+e*cAP)*ci,
                    sAN*(sth+e*sAP) - cAN*(cth+e*cAP)*ci,
                    -(cth+e*cAP)*si};
    ];
    
    Return[{rVec,vVec}];
]



rv2elem[mu_, rVec_, vVec_] := Block[{a,e,i,AN,AP,f,eps,ir,ie,ip,ih},
    
    (* define a small number *)
    eps = 0.000000000001;
    
    (* compute orbit radius *)
    r = Norm[rVec];
    ir = rVec/r;
    
    (* compute the angular momentum vector *)
    hVec = Cross[rVec,vVec];
    h    = Norm[hVec];
    
    (* compute the eccentricity vector *)
    cVec = Cross[vVec, hVec] - mu*rVec/r;
    e    = Norm[cVec]/mu;
    
    (* compute semi-major axis *)
    ai = 2./r - (Dot[vVec,vVec])/mu;
    If[ (Abs[ai] > eps),
        (* elliptic or hyperbolic case *)
        a = 1/ai;
    ,
        (* parabolic case *)
        p  = h^2/mu;
        rp = p/2;
        a  = -rp;   (* a is not defined for parabola, so -rp is returned instead *)
        e  = 1;
    ];
    
    If[ (h < eps) ,  (* rectilinear motion case *)
        ie = ir;
        (* ip and ih are arbitrary *)
        dum  = Cross[ie,{0,0,1}];
        dum2 = Cross[ie,{0,1,0}];
        If[ (Norm[dum] > Norm[dum2]) , 
            ih = dum/Norm[dum];
        ,
            ih = dum2/Norm[dum2];
        ];
        ip = Cross[ih, ie];
    ,
        (* compute perifocal frame unit direction vectors *)
        ih = hVec/h;
        If[ (Abs[e] > eps),
            (* non-circular case *)
            ie = cVec/mu/e;
        ,
            (* circular orbit case.  Here ie, ip are arbitrary, as long as they
               are perpenticular to the ih vector.  *)
            ie = ir;
        ];
        ip = Cross[ih,ie];
    ];
    
    (* compute the 3-1-3 orbit plane orientation angles *)
    AN = ArcTan[-ih[[2]],ih[[1]]];
    i  = ArcCos[ih[[3]]];
    AP = ArcTan[ip[[3]],ie[[3]]];
    
    
    If[ (h < eps) ,          (* rectilinear motion case *)
        If[ (ai > 0) ,       (* elliptic case *)
            Ecc = ArcCos[1-r*ai];
            If[ (Dot[rVec,vVec] > 0) ,
                Ecc = 2*\[Pi] - Ecc;];
            f = Ecc;    (* for this mode the eccentric anomaly is returned *)
        ,            
            (* hyperbolic case *)
            H = ArcCosh[(r*ai+1)];
            If[ (Dot[rVec,vVec] < 0) ,
                H = 2*\[Pi] - H;];
            f = H;    (* for this mode the hyperbolic anomaly is returned *)
        ];
    ,
        (* compute true anomaly *)
        f = ArcTan[Dot[ie,ir],Dot[Cross[ie,ir],ih]];
    ];

    
    Return[{a,e,i,AN,AP,f}];
]



AtmosphericDensity[alt_] := Block[{logdensity, density, val},
    
    (* Smooth exponential drop-off after 1000 km *)
    If [alt > 1000.,
        logdensity = (-7*10^(-05))*alt-14.464;
        density = 10^logdensity;
        Return[density];
    ];
     
    (* Calculating the density based on a scaled 6th order polynomial fit 
       to the log of density *)
    val = (alt-526.8000)/292.8563;
    logdensity = 0.34047*val^6 -0.5889*val^5 -0.5269*val^4 + 1.0036*val^3 + 0.60713*val^2 -2.3024*val -12.575;
    
    (* Calculating density by raising 10 to the log of density *)
    density = 10^logdensity;
    
    Return[density];
    
]





AtmosphericDrag::"altitudeCheck" = " AtmosphericDrag() received rvec = `1`.  The value of rvec should produce a positive altitude for the Earth."
AtmosphericDrag[Cd_, A_, m_, rvec_, vvec_] := Block[{alt, r, v, ad, density},

    (* find the altitude and velocity *)
    r   = Norm[rvec];
    v   = Norm[vvec];
    alt = r - astroConstants`reqEarth;
    
    
    (* Checking if user supplied a orbital position is inside the earth *)
    If [(alt <= 0.) ,
        advec= { Null, Null, Null};
        Message[AtmosphericDrag::"altitudeCheck", rvec];
        Return[advec];
    ];
    
    
    (* get the Atmospheric density at the given altitude in kg/m^3 *)
    density = AtmosphericDensity[alt];
    
    (* compute the magnitude of the drag acceleration *)
    ad=((-0.5)*density*(Cd*A/m)*((v*1000)^2))/1000.;
    
    (* computing the vector for drag acceleration *)
    advec = ad/v*vvec;
    
    
    Return[advec];
]


JPerturb::"OrderCheck" = " JPerturb() received order = `1`.  The order must be a value between [2,6]."
JPerturb[rvec_, num_] := Block[{x,y,z,mu,req,r},
    
    (* Constants for Earth *)
    mu  = astroConstants`\[Mu]Earth;
    req = astroConstants`reqEarth;
    J2 = astroConstants`J2;
    J3 = astroConstants`J3;
    J4 = astroConstants`J4;
    J5 = astroConstants`J5;
    J6 = astroConstants`J6;
        
    (* Calculate the J perturbations *)
    x=rvec[[1]];
    y=rvec[[2]];
    z=rvec[[3]];
    r=Norm[rvec];
    
    (* Error Checking *)
    If[ ((num < 2) || (num > 6)) ,
        Message[JPerturb::"OrderCheck", num];
        ajtot = {Null,Null,Null};
        Return;
    ];
    
    (* Calculating the total acceleration based on user input *)
    If[ num >= 2,
        ajtot=(-3/2*J2*(mu/(r^2))*(req/r)^2)*{(1-5*(z/r)^2)*(x/r),(1-5*(z/r)^2)*(y/r),(3-5*(z/r)^2)*(z/r)};
    ];
    If[ num >= 3,
        ajtot=ajtot+(1/2*J3*(mu/(r^2))*(req/r)^3)*{5*(7*(z/r)^3-3*(z/r))*(x/r),5*(7*(z/r)^3-3*(z/r))*(y/r),3*(-10*(z/r)^2+(35/3)*(z/r)^4+1)};
    ];
    If[ num >= 4,
        ajtot=ajtot+(5/8*J4*(mu/(r^2))*(req/r)^4)*{(3-42*(z/r)^2+63*(z/r)^4)*(x/r),(3-42*(z/r)^2+63*(z/r)^4)*(y/r),(15-70*(z/r)^2+63*(z/r)^4)*(z/r)};
    ];
    If[ num >= 5,
        ajtot=ajtot+(1/8*J5*(mu/(r^2))*(req/r)^5)*{3*(35*(z/r)-210*(z/r)^3+231*(z/r)^5)*(x/r),3*(35*(z/r)-210*(z/r)^3+231*(z/r)^5)*(y/r),-(15-315*(z/r)^2+945*(z/r)^4-693*(z/r)^6)};
    ];
    If[ num >= 6,
        ajtot=ajtot+(-1/16*J6*(mu/(r^2))*(req/r)^6)*{(35-945*(z/r)^2+3465*(z/r)^4-3003*(z/r)^6)*(x/r),(35-945*(z/r)^2+3465*(z/r)^4-3003*(z/r)^6)*(y/r),-(3003*(z/r)^6-4851*(z/r)^4+2205*(z/r)^2-245)*(z/r)};
    ];

    Return[ajtot];
]






End[]

EndPackage[]
