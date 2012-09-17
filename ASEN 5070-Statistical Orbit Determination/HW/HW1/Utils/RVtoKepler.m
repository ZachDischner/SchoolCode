% Function to convert cartesian coordinates to Keplarian
% 
% Written by Zach Dischner on 8/28/2012
% 
% inputs:
%         R: [i j k] radius vector     km
%         V: [i j k] r' vector         km/s
%         u: gravatational constant    km^3/s^2
%         
% outputs:
%         a: Semi-Major Axis           km
%         e: Eccentricity              []
%         i: Inclination               Rad
%         Omega: Right Ascension       Rad
%         w: 
%         nu: True Anomoly             Rad
%         EA: Eccentric Anomoly        Rad
%         T: Time of Periapse passage
%         
%         
%        

function [a e i Omega w nu EA T] = RVtoKepler(R,V,u,t)

% If no provided time input
if ~exist('t','var')
    t = 0;
end

%% 1.0 Compute Specific Angular Momentum
H    = cross(R,V);      %units?
h    = norm(H);

%% 2.0 Compute Scalar Radius and Velocity
r    = norm(R);             % Km
v    = norm(V);             % Km/s

%% 3.0 Compute Specific Energy E, verify Elliptical Motion
E    = (v).^2./2 - u./r;    % units?
% Verify elliptical? Look it up. 

%% 4.0 Compute Semimajor Axis, a
a   = -u./(2*E);            % km

%% 5.0 Compute Eccentricity, e
e   = sqrt(1-((h.^2)./(a.*u)));   % []

%% 6.0 Compute Inclination,             0 < i < 180
i   = acos(H(3)./h);             % Rad

%% 7.0 Compute Right Ascension of the Ascending Node,       0 < Omega < 360
Omega   = atan2(H(1),-H(2)); % Rad

%% 8.0  True Anomoly           0 < nu < 360
p   = a.*(1-e.^2);
nu  = atan2((sqrt(p./u)).*(dot(V,R)),p-r); 

%% 9.0 Argument of Periapse    0 < w < 360
w = mod(atan2(R(3)./sin(i),(R(1).*cos(Omega) + R(2).*sin(Omega))) - nu,pi);

%% 10.0 Compute Eccentric Anomoly,   0 < EA < 360
EA = 2*atan( sqrt((1-e)./(1+e)).*tan(nu./2));  % Radians
EA = unwrap(EA);

%% 11.0 Compute Time of Periapse Passage, T
n = sqrt(u./(a.^3));
T = t - (1./n) .* (EA - e.*sin(EA));



%% Handle all Conversions (possible add as parser input)
% EA = EA.*180./pi;











end