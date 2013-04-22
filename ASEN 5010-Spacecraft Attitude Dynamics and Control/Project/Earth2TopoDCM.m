function TE = Earth2TopoDCM( lambda, phi)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           Earth2TopoDCM.m
% Author:   Zach Dischner
% Date:     April 4, 2013
% 
% Usage:
%   TE = Earth2TopoDCM( lambda, phi, inclination )
%
% Description:  Maps the earth-fixed frame [EN] to a topographic frame 
%               [TE] ==> {n e d}
% 
%               Lamens ==> Rotation matrix between earth-fixed and
%                          topographic frame
% 
% Inputs:  lambda ==> Angle between grenwich and equator orbit crossing
%          phi    ==> Along-track angle between satellite and equator plane
%                     (satellite latitude)
% ???      inclination ==> Orbit inclination?
%
% Outputs: TE     ==> DCM from earth fixed frame to the topographic frame
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%% Define Individual Rotations
% First: Rotate LAMBDA about e^3
rot1    = Euler3(lambda);

% Second: Rotate (-PHI  - pi/2)about e^2' (Can be done in just one
% rotation.
rot2    = Euler2(-phi - pi/2);

% Third: Inclination angle? 
inclination = 0;
rot3    = Euler1( inclination );

%% Assemble the Rotation Matrix
TE = rot3*rot2*rot1;

