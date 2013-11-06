%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        calculateSatellitePosition.m
% Author      : Zach Dischner
% Date        : 10/24/2013
% Description : calculate satellite position from GPS ephemeris dataset
% 
%
%                 __...____________________          ,
%                `(\ [ ===NCC-1700===--|__|) ___..--"_`--.._____
%                  `"""""""""""""""""| |""` [_""_-___________"_/
%                                    | |   /..../`'-._.-'`
%                                ____| |__/::..'_
%                               |\ ".`"` '_____//\
%                               `"'-.   """""  \\/
%                                    `""""""""""`
% Inputs        : ephem - Satellite ephemeris dataset
%                 params - structure containing keplarian specs and extra
%                           calculations for sat position
% 
% Outputs       : [rk]-3d ECI coordinates of satellite
% 
% History       October 11 2013 - First Rev
%               October 24 2013 - Reformatted output to [rk,tk]
%               October 30 2013 - Added params to return, check for a 
%                                   better Ek guess (to speed up 'find')
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [rk,dt_rel, params] = calculateSatellitePosition(ephem,params)

%-----Extract all ephemeris components to make life easy and epicer
ephem = num2cell(ephem);
[prn,M0,delta_n,ecc,sqrt_a,Loa,incl,perigee,ra_rate,i_rate,Cuc,Cus,Crc,Crs,Cic,Cis,...
    Toc,IODE,GPS_week,Toc,Af0,Af1,Af2,nil,health] = deal(ephem{:});
A = sqrt_a^2;
%------Correct Mean Motion
n0  = sqrt(params.mu/(A)^3); % Calculated mean motion [rad/s]
n   = n0 + delta_n;                 % Corrected Mean Motion

%------Correct Time
tk  = params.Secs-Toc;  

%------Mean Anomaly
Mk = M0 + n*tk; % Mean anomaly

%------Eccentric Anomaly
if isfield(params,'E_guess')
    guess=params.E_guess;
else
    guess=0;
end
Ek = fsolve(@(Ek) (Ek)-ecc*sin(Ek)-Mk,guess,params.options);
params.E_guess = Ek;

%------True Anomaly
vk = atan2(      (sqrt(1-ecc^2)*sin(Ek)/(1-ecc*cos(Ek))), ...
                      ((cos(Ek)-ecc)/(1-ecc*cos(Ek)))   );
                  
%------Argument of Latitude
Phik = vk + perigee;

%------Second Harmonic Perturbations
del_uk = Cus*sin(2*Phik) + Cuc*cos(2*Phik);
del_rk = Crs*sin(2*Phik) + Crc*cos(2*Phik);
del_ik = Cis*sin(2*Phik) + Cic*cos(2*Phik);

%------Corrected argumet of latitude, radius, inclination
uk = Phik + del_uk;
rk = A*(1-ecc*cos(Ek)) + del_rk;
ik = incl + del_ik + i_rate*tk;

%------Position in Orbit Plane
xkp = rk*cos(uk);
ykp = rk*sin(uk);

%------Corrected Longitude of ascending node
Omegak = Loa + (ra_rate - params.we)*tk - params.we*Toc;

%------Earth Fixed Coordinates
xk = xkp * cos(Omegak) - ykp * cos(ik) * sin(Omegak);
yk = xkp * sin(Omegak) + ykp * cos(ik) * cos(Omegak);
zk = ykp * sin(ik);

%------Relativity time shift
dt_rel = 2*sqrt(params.mu)/params.c^2 * ecc * sqrt_a * sin(Ek);
rk = [xk,yk,zk];

















