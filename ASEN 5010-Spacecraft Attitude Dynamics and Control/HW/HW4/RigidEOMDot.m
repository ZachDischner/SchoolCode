%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           RigidEOMDot.m
% Author:   Zach Dischner
% Date:     March 19, 2013
% 
% Usage:
%   wdot = RigidEOMDot(w)
%
% Description:  Used for integrating EOM equations for a rigid body by 
%               determining the derivative of angular velocity
% 
% Inputs:  X    => [Angular velocity vector;321Vector]
%          I    => Principle Inertia Matrix
%
% Outputs: Xdot => Derivative of angular velocity vector
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

function Xdot = RigidEOMDot(t,X,I)

w   = X(1:3);
ang = X(4:6); 
I = diag(I);

wdot(1) = -(I(3) - I(2))*w(2)*w(3)/I(1);
wdot(2) = -(I(1) - I(3))*w(3)*w(1)/I(2);
wdot(3) = -(I(2) - I(1))*w(1)*w(2)/I(3);


angdot   = BmatEuler321(ang)*w';

Xdot=[wdot';angdot];



end