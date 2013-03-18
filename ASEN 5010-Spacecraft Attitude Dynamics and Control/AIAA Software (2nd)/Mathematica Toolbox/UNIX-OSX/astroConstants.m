(* :Title: astroConstants *)

(* :Author: Hanspeter Schaub *)

(* :Summary:
This package provides common constants of astrydynamics problems.

	
*)


(* :Package Version: 1.0. *)

(* :Copyright: Freeware *)

(* :History: 
	Version 1.0 by Hanspeter Schaub (Virginia Tech), June 2005.
*)

(* :Keywords: orbital motion, astrodynamics *)

(* :Mathematica Version: 5.0 *)

BeginPackage["astroConstants`"]

(* 
  universial Gravitational Constant
  units are in km^3/s^2/kg
*)
Guniversial = 6.67259e-20       ;   

(* 
   astronomical unit in units of kilometers
*)
AU          = 149597870.691     ;


(* 
    common conversions
*)
D2R         = \[Pi]/180.;
R2D         = 180./\[Pi];

(*
  Gravitational Constants mu = G*m, where m is the planet of the attracting
  planet.  All units are km^3/s^2. 
  values are obtained from http://ssd.jpl.nasa.gov/astro_constants.html
*)
\[Mu]Sun       = 1.32712440018*10^(11)  ;
\[Mu]Mercury   = 2.2032*10^(4)          ;
\[Mu]Venus     = 3.2485859*10^(5)       ;
\[Mu]Earth     = 3.986004418*10^(5)     ;
\[Mu]Moon      = 4.9027988*10^(3)       ;
\[Mu]Mars      = 4.28283*10^(4)         ;
\[Mu]Jupiter   = 1.2671277*10^(8)       ;
\[Mu]Saturn    = 3.79406*10^(7)         ;
\[Mu]Uranus    = 5.79455*10^(6)         ;
\[Mu]Neptune   = 6.83653*10^(6)         ;
\[Mu]Pluto     = 983.              ;



(*
 planet information for major solar system bodies. Units are in km.
  data taken from http://nssdc.gsfc.nasa.gov/planetary/planets.html
*)
(* Sun: *)
reqSun      = 695000.           ;

(* Mercury:  *)
reqMercury  =   2439.7          ;
J2Mercury   =     60.0*10^(-6)  ;
smaMercury  =      0.38709893*AU;
iMercury    =      7.00487*D2R  ;
eMercury    =      0.20563069   ;

(* Venus: *)
reqVenus    =   6051.8          ;
J2Venus     =      4.458*10^(-6);
smaVenus    =      0.72333199*AU;
iVenus      =      3.39471*D2R  ;
eVenus      =      0.00677323   ;

(* Earth:*)
reqEarth    =   6378.14         ;
smaEarth    =      1.00000011*AU;
iEarth      =      0.00005*D2R  ;
eEarth      =      0.01671022   ;

(* Moon:*)
reqMoon     =   1738.1          ;
J2Moon      =    202.7*10^(-6)  ;
smaMoon     =      0.3844e6     ;
eMoon       =      0.0549       ;

(* Mars:*)
reqMars     =   3397.2          ;
J2Mars      =   1960.45*10^(-6) ;
smaMars     =      1.52366231*AU;
iMars       =      1.85061*D2R  ;
eMars       =      0.09341233   ;

(* Jupiter:*)
reqJupiter  =  71492.           ;
J2Jupiter   =  14736.*10^(-6)   ;
smaJupiter  =      5.20336301*AU;
iJupiter    =      1.30530*D2R  ;
eJupiter    =      0.04839266   ;

(* Saturn:*)
reqSaturn   =  60268.           ;
J2Saturn    =  16298.*10^(-6)   ;
smaSaturn   =      9.53707032*AU;
iSaturn     =      2.48446*D2R  ;
eSaturn     =      0.05415060   ;

(* Uranus: *)
reqUranus   =  25559.           ;
J2Uranus    =   3343.43*10^(-6) ;
smaUranus   =     19.19126393*AU;
iUranus     =      0.76986*D2R  ;
eUranus     =      0.04716771   ;

(* Neptune: *)
reqNeptune  =  24746.           ;
J2Neptune   =   3411.*10^(-6)   ;
smaNeptune  =     30.06896348*AU;
iNeptune    =      1.76917*D2R  ;
eNeptune    =      0.00858587   ;

(* Pluto: *)
reqPluto    =   1137.           ;
smaPluto    =     39.48168677*AU;
iPluto      =     17.14175*D2R  ;
ePluto      =      0.24880766   ;    


(*
  Zonal gravitational harmonics Ji for the Earth
*)
J2 = 1082.63*10^(-6);
J3 =   -2.52*10^(-6);
J4 =   -1.61*10^(-6);
J5 =   -0.15*10^(-6);
J6 =    0.57*10^(-6);







EndPackage[]
