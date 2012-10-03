function [yprime] = SpringDeriv(time,y)

% Simple function to pass to integrator for HW5, problem 3  => 4.15 in book. 
% yprime = d(y)/dt = 

% where y = [x v PHI]
% 
% so yprime = [    0    1       [x
%                -w^2   0]       v]     ; A*Phi
% 
% 
% where w^2 = (k1 + k2)/m

k1          = 2.5;  % N/m
k2          = 3.7;  % N/m
m           = 1.5;  % Kg
w_squared   = (k1 + k2)/m;

A           = [ 0 1; -w_squared 0];

yprime =[ [0 1; -w_squared 0]*y(1:2); reshape( A * reshape(y(3:6),2,2), 4,1)];






end
