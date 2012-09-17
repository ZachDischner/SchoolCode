% Function to convert cartesian coordinates to Keplarian
% 
% Written by Zach Dischner on 8/28/2012
% 
% inputs:
%         a: Semi-Major Axis           km
%         e: Eccentricity              []
%         i: Inclination               Degrees
%         Omega: Right Ascension       Degrees
%         w: Argument of Periapses     Degrees 
%         nu: True Anomoly             Degrees
%         u: gravatational constant    km^3/s^2
%        (t): optional, time           s
%      
%         
% outputs: (All in ECI I think)
%         R: [i j k] radius vector     km
%         V: [i j k] r' vector         km/s
%         
%         
%

function[R, V] = Kepler_to_RV_easy(a, e, i, Omega, w, nu, u)

p = a*(1 - e^2);                           % Semiparameter (Km)

%%  Perifocal Position and Velocity
R_per = [(p.*cos(nu)) / (1 + e.*cos(nu));    % X (Km)
         (p.*sin(nu)) / (1 + e.*cos(nu));    % Y (Km)
                      0;                     % Z (Km)
        ];
 
V_per = [ -(u./p).^(1/2) .* sin(nu);         % X' (Km/s)
         (u./p).^(1/2) .* (e + cos(nu));     % Y' (Km/s)
                      0;                     % Z' (km/s)
        ];
                    
%%  Transformation Matrix (3 Rotations)  %%%
Per2ECI = [ cos(Omega)*cos(w)-sin(Omega)*sin(w)*cos(i) ...
           (-1)*cos(Omega)*sin(w)-sin(Omega)*cos(w)*cos(i) ...
                       sin(Omega)*sin(i); ...
             sin(Omega)*cos(w)+cos(Omega)*sin(w)*cos(i) ...
             (-1)*sin(Omega)*sin(w)+cos(Omega)*cos(w)*cos(i) ...
                       (-1)*cos(Omega)*sin(i); ...
               sin(w)*sin(i)  cos(w)*sin(i)  cos(i)];

%% Final Step, Perifocal to ECI Coordinate Transformation 
R = Per2ECI*R_per;
V = Per2ECI*V_per;








