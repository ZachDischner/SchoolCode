function [yprime] = RV_Deriv1a(time,y)

% Simple function to pass to integrator for HW5, problem 1. 
% yprime = d(y)/dt = 

% where y = [x y x' y']
% 
% so yprime = [x' y' x'' y'']
% 
% where x'' = -x/r^3;  y'' = -y/r^3;

u = 1;

r = sqrt(y(1).^2 + y(2).^2);
      
yprime = [y(3) y(4) -y(1)/r^3 -y(2)/r^3 ]';

end