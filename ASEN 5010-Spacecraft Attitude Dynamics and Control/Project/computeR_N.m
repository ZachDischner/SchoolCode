function r_N = computeR_N(t)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           computeR_N.m
% Author:   Zach Dischner
% Date:     April 1, 2013
% 
% Usage:
%   r_N = computeR_N(t)
%
% Description:  Computes orbit position vector in inertial frame components
%               for ASEN 5010 project circular orbit
% 
% Inputs:  t    ==> Time value array [s]. Either row or column vector
%
% Outputs: r_N  ==> Orbit position in N frame components
%                   | r1(t1)    r2(t1)  r3(t1) |
%                   | r1(t2)    r2(t2)  r3(t2) |
%                   | r1(t3)    r2(t3)  r3(t3) |
%                                ...
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Orbit Parameters
Omega   = deg2rad(20);      % RAAN?                     [rad]
i       = deg2rad(75);      % Orbit Inclination Angle   [rad]
mu      = 398600;           % Gravatational Parameter   [km^3/s^2]
r       = 6878;             % Circular Orbit Radius     [km]
n       = sqrt(mu/r^3);     % Mean Orbit Rate           [s]
theta_0 = 0;                % Initial Angle             [rad]

%% Compute r_N(t)
if ~iscolumn(t)
    t=t';
end

theta   = theta_0 + n*t;

r_N = r*[ cos(Omega).*cos(theta)-sin(Omega).*sin(theta).*cos(i) , ...
        sin(Omega).*cos(theta)+cos(Omega).*sin(theta).*cos(i) , ...
                            sin(theta).*sin(i)                  ];
        

