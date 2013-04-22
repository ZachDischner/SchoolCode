function EN = ECI2ECF( t )
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           ECI2ECF.m
% Author:   Zach Dischner
% Date:     April 6, 2013
% 
% Usage:
%   EN = ECI2ECF( t )
%
% Description:  Maps the earth-centered inerital frame [N] to an
%               earth-centered earth-fixed frame [E]. Assumes constant
%               rotational rate
% 
%               Lamens ==> Rotates ECI to ECF by [EN]
% 
% Inputs:  t  ==> Time that has passed since grenwich 0:00:00 ???
%
% Outputs: EN ==> DCM from ECI to ECEF.
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Obtain Rotation Rate
w_EN = 361 *(pi/180)*(1/24)*(1/3600); %rotation [deg/day][rad/deg][day/hour][hour/sec]

%% Obtain Rotation Angle
gamma_t0 = deg2rad(20);

gamma    = gamma_t0 + w_EN*t;

%% Assemble DCM
EN       = Euler3(gamma);

