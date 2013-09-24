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

function Xprime = diff321w_part4(t,X,I,P,K)

% w   = X(4:6)';
% sigma = X(1:3)'; 
% 
% 
% wdot        = inv(I) *( - tilde(w)*I*w-P*w-K*sigma );
% 
% sigmadot    = 1/4*BmatMRP( sigma )*w;
% 
% Xdot=[sigmadot;wdot];

Xprime = [   zeros(3,3)    ,    BmatMRP(X(1:3)')*1/4     ; ...
            -inv(I)*K  ,   -inv(I)*P]*X';


end