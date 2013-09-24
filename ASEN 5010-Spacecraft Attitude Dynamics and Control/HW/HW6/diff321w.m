%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           diff321W.m
% Author:   Zach Dischner
% Date:     March 19, 2013
% 
% Usage:
%   Xdot = diff321w(t,X, I)
%
% Description:  Used for integrating EOM equations for a rigid body by 
%               determining the derivative of angular velocity
% 
% Inputs:  X    => [Angular velocity vector;321Vector]
%          I    => Principle Inertia Matrix
%          U    => Handle for control function ( U =@(x,y,z) ...)
%          L    => External torque
%
% Outputs: Xdot => Derivative of angular velocity vector
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

function Xdot = diff321w(t,X,I,L,U,K,P)

w   = X(1:3)';
ang = X(4:6)'; 



wdot = inv(I) *( - tilde(w)*I*w + U(ang,w,L,K,P) + L);

angdot   = BmatEuler321(ang)*w;

Xdot=[wdot;angdot];

end