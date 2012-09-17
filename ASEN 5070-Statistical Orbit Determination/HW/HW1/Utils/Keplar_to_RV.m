% Function to convert cartesian coordinates to Keplarian
% 
% Written by Zach Dischner on 8/28/2012
% 
% inputs:
%         a: Semi-Major Axis           km
%         e: Eccentricity              []
%         i: Inclination               Degrees
%         Omega: Right Ascension       Degrees
%         w: Argument of Periapses     Degrees 
%         nu: True Anomoly             Degrees
%         T: Time since Periapses
%         u: gravatational constant    km^3/s^2
%        (t): optional, time           s
%      
%         
% outputs:
%         R: [i j k] radius vector     km
%         V: [i j k] r' vector         km/s
%         
%         
%        

function [R V] = RVtoKepler(a,e,i,Omega,w,T,u,t)

if ~exist('t','var')
    t = 0;
end

%% 1.0 Compute Mean Anomoly
n = sqrt(u./(a.^3));
M = n.*(t - T);

%% 2.0 compute Eccentric Anomoly EA       0 < EA < 360
EA = fsolve( @(EA) EA - e.*sin(EA) - M,pi);    % Radians 

%% 3.0 Compute True Anomoly               0 < nu < 360
nu = 2*atan2( ( sqrt( (1 + e)./(1 -e)) .* tan(EA/2)));

nu = nu.*180/pi;

%% 4.0 Compute Radius, r
r = a.*(1-e.*cos(EA));    % Or   r = a(1-e^2)/(1+e*cos(nu))

%% 5.0 Compute Specific Angular Momentum
h = sqrt( u.*a.*( 1-e.^2 ) );

%% 6.0 Compute Position components   R => [Xi Yj Zk]
R = [r.*( cosd(Omega).*cosd(w + nu) - sind(Omega).*sind(w + nu).*cosd(i));...
     r.*( sind(Omega).*cosd(w + nu) + cosd(Omega).*sind(w + nu).*cosd(i));...
     r.*( sind(i).*sind(w + nu))...
     ];
 
%% 7.0 Compute Velocity Components   V => [X'i Y'j Z'k]






end