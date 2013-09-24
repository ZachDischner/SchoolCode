%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           diffMRPw.m
% Author:   Zach Dischner
% Date:     March 19, 2013
% 
% Usage:
%
% Description:  Used for integrating EOM equations for a rigid body by 
%               determining the derivative of angular velocity
% 
% Inputs:  X    => [Angular velocity vector;MRP]
%          I    => Principle Inertia Matrix
%          U    => Handle for control function ( U =@(x,y,z) ...)
%          L    => External torque
%          K/P  => Gain matrix for U
%
% Outputs: Xdot => Derivative of angular velocity vector
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

function Xdot = diff321w(t,X,I,L,U,K,P,K1)

w   = X(1:3)';
sigma = X(4:6)'; 

z = X(7:9);

wdot        = inv(I) *( - tilde(w)*I*w + U(sigma,w,L,K,P,K1,z') + L);

sigmadot    = 1/4*BmatMRP( sigma )*w;

Xdot=[wdot;sigmadot; sigma];

end