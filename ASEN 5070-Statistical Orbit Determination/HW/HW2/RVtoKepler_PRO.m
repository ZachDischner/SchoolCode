% Name:      RVtoKepler.m
% Author:    Jeffrey S. Parker
% Purpose:   A function used to convert an object's 
%            Cartesian coordinates to Keplerian elements
%
% Call:      [a, E, i, Omega, w, nu] = RVtoKepler(Rxyz, Vxyz, mu)
%
% Inputs:    Rxyz         = (x,y,z) Cartesian Position (km)
%            Vxyz         = (vx,vy,vz) Cartesian Velocity (km/s)
%            mu           = Gravitational parameter of the system (GM), km^3/s^2
%
% Outputs:   a            = Semi-major axis (km)
%            E            = Eccentricity
%            i            = Inclination (rad)
%            Omega        = Right-Ascension of the ascending node (rad)
%            w            = Argument of periapse (rad)
%            nu           = True anomaly (rad)
%
% Required Subroutines:   N/A

function[a, E, i, Omega, w, nu] = RVtoKepler_PRO(Rxyz, Vxyz, mu)

R = norm(Rxyz,2)           % Magnitude of R
V = norm(Vxyz,2)           % Magnitude of V

% Compute the Specific Angular Momentum Vector
h = cross(Rxyz,Vxyz) 
H = norm(h,2)              % Magnitude of h
if (H == 0)
    fprintf('Warning!  Specific Angular Momentum (H) = 0 \n')
    warning off;
end
k = [0 0 1];                % k-unit vector

n = cross(k,h)             % line of nodes
N = norm(n,2)              % magnitude of line

Energy = (1/2)*V^2 - mu/R  % Specific Energy

a = (-1/2)*mu / Energy;     % semi-major axis

e = cross(Vxyz,h)./mu - Rxyz./R; % eccentricity vector
E = norm(e,2);              % eccentricity
% or we can use: E = sqrt(1+2*Energy*H^2/mu^2)

i = acos(h(3)/H);           % inclination
while (i < 0)               % Quadrant checks
    i = i + pi;
end
while (i >= pi)
    i = i - pi;
end

%%%  Special Cases Check  %%%
if (i < 0.005 & E >= 0.005)
    fprintf('\n Warning!  Inclination very close, or equal, to zero! \n')
    fprintf('(Equatorial Orbit) \n\n')
    warning off;
    w_true = acos(e(1)/E);  % True Longitude of Periapse
    if (e(2) < 0)           % Quadrant Check
        w_true = 2*pi - w_true;
    end
    fprintf('If you would like it, the True Longitude of Periapse \n')
    fprintf('(omega_true) = %.8f \n',w_true)
end

if (E < 0.005 & i > 0.005)
    fprintf('\n Warning!  Eccentricity very close, or equal, to zero! \n')
    fprintf('(Circular Orbit) \n\n')
    warning off;
    u = acos(dot(n,Rxyz) / (N*R));  % Argument of Latitude
    if (Rxyz(3) < 0)                % Quadrant Check
        u = 2*pi - u;
    end
    fprintf('If you would like it, the Argument of Latitude (u) = %.8f \n',u)
end

if (i < 0.005 & E < 0.005)
    fprintf('\n Warning!  Both the inclination and eccentricity \n')
    fprintf('are very close, or equal, to zero! \n')
    fprintf('(Circular, Equatorial Orbit) \n\n')
    lambda = acos(Rxyz(1)/R);   % True Longitude
    if (Rxyz(2) < 0)            % Quadrant Check
        lambda = 2*pi - lambda;
    end
    fprintf('If you would like it, the True Longitude \n')
    fprintf('(lambda_true) = %.8f \n',lambda)
end


Omega = acos(n(1)/N)       % Right Ascension of the Ascending Node
if (n(2) < 0)               % Quadrant Check
    Omega = 2*pi - Omega
end

w = acos(dot(n,e) / (N*E)); % Argument of Periapse
if (e(3) < 0)               % Quadrant Check
    w = 2*pi - w;
end

nu = acos(dot(e,Rxyz) / (E*R)); % True Anomaly
if (dot(Rxyz,Vxyz) < 0)         % Quadrant Check
    nu = 2*pi - nu;
end

warning on;