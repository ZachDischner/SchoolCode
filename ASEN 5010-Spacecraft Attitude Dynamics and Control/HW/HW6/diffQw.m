%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           diffQw.m
% Author:   Zach Dischner
% Date:     March 19, 2013
% 
% Usage:
%
% Description:  Used for integrating EOM equations for a rigid body by 
%               determining the derivative of angular velocity
% 
% Inputs:  X    => [Angular velocity vector;B]
%          I    => Principle Inertia Matrix
%          U    => Handle for control function ( U =@(x,y,z) ...)
%          L    => External torque
%          K/P  => Gain matrix for U
%
% Outputs: Xdot => Derivative of angular velocity vector
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

function Xdot = diff321w(t,X,I,L,U,K,P)

w   = X(1:3)';
B = X(4:7)'; 

epsilon = B(2:4);
wdot        = inv(I)*( - tilde(w)*I*w + U(epsilon,w,L,K,P) + L);

sigmadot    = 1/2*BmatEP( B )*w;

Xdot=[wdot;sigmadot];

end