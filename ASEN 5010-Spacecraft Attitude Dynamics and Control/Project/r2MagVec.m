function M_N = r2MagVec( r_N, t )
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           r2MagVec.m
% Author:   Zach Dischner
% Date:     April 6, 2013
% 
% Usage:
%   M_N = r2MagVec( r_N, t )
%
% Description:  Compute the Magnetic field vector for a given satellite
%               position and time
% 
%               Lamens ==> Find the magnetic field vector at a point in
%                          space and time
% 
% Inputs:   r_ECEF  ==> ECEF satellite position vector [r1 r2 r3]
%           t    ==> Time that has passed since grenwich 0:00:00 ???
%
% Outputs:  M_N  ==> Magnetic field vector in N frame components
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Defines and Constants
r_E     = 6378; % Earth's Equatorial Radius [km]
r_orbit = norm(r_N);
% lam     = 0.001*sin(t);
% phi     = 0.1*cos(t);
lam     = atan2(r_N(2), r_N(1)); % *look at it again
phi     = asin(r_N(3)/r_orbit);

%% Calculate Magnetic Field Vector
M_T     = -(r_E/r_orbit)^3.*...
                        [ -cos(phi)     sin(phi)*cos(lam)     sin(phi)*sin(lam);...
                            0             sin(lam)              -cos(lam)   ;...
                       -2*sin(phi)  -2*cos(phi)*cos(lam)   -2*cos(phi)*sin(lam)]...
                       *[29900; 1900; -5530];       % nT
                   
                   
EN      = ECI2ECF(t);
TE      = Earth2TopoDCM(lam,phi);

M_N     = EN'*TE'*M_T;

M_N = M_N/norm(M_N);

