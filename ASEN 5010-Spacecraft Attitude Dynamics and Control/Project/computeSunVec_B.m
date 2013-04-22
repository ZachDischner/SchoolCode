function s_B = computeSunVec_B(sigma_BN)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           computeSunVec_B.m
% Author:   Zach Dischner
% Date:     April 4, 2013
% 
% Usage:
%   s_B = computeSunVec_B(sigma_BN)
%
% Description:  Computes the sun attitude attitude as seen by the satellite 
%               body. It does so with a constant inertial sun attitude, and
%               the MRP set sigma_B/N.
% 
%               Given s_N, and BN, get s_B
% 
% Inputs:  sigma_BN    ==> MRP set describing rotation between B and N
%                          frames
%
% Outputs: s_B  ==> Sun direction vector in B frame components
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Data Masseuse 
if ~iscolumn(sigma_BN)
    sigma_BN = sigma_BN';
end

%% Define Sun Position Vector
s_N = [ 0 ; -1 ; 0 ];

%% Convert MRPs into [BN] DCM
BN  = MRP2C(sigma_BN);

%% Get the Rotation
s_B = BN*s_N;

s_B = s_B/norm(s_B);
