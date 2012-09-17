% MatRV2Keplar
% 
% Function to convert cartesian coordinates to Keplarian elements, if the
% R,V vectors are matrices. 
% 
% Written by Zach Dischner on 9/6/2012, with help from Greg Nelson for
% helping debugging 
% 
% inputs:
%         R: [i j k] radius vector     km
%         V: [i j k] r' vector         km/s
%         u: gravatational constant    km^3/s^2
%         t: (optional) time 
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

function [a e i Omega w v tp]=MatRV2Kepler(R,V,u,t)

% Is time defined?
if ~exist('t','var')
    t = 0;
end
[m n]=size(R);

%% Specific Angular Momentum
h = cross(R,V);

%% R,V,H magnitudes
h_mag = zeros(length(m));
R_mag = h_mag;
V_mag = h_mag;
for ii=1:m
    h_mag(ii,1) = norm(h(ii,:));
    R_mag(ii,1) = norm(R(ii,:));
    V_mag(ii,1) = norm(V(ii,:));
end

%% Specific Energy
E = V_mag.^2/2-u./R_mag;

%% Semi-major axis
a = -u./(2*E);

%% Eccentricity
e = sqrt(1-h_mag.^2./(a*u));

%% Inclination
i = acos(h(:,3)./h_mag);

%% RAAN
Omega = atan2(h(:,1),-h(:,2));

%% Argument of Periapsis 
% semi-parameter
p = a.*(1-e.^2);

% Calculate true anomaly
for ii=1:m
    v(ii,1) = atan2(sqrt(p(ii)/u)*dot(V(ii,:),R(ii,:)),p(ii)-R_mag(ii));
end

% Argument of periapsis 
w = mod(atan2(R(:,3)./sin(i),R(:,1).*cos(Omega)+R(:,2).*sin(Omega)) - v,pi);

% Calculate argument of periapsis
% w = argLat-v;
% w = mod(w,pi);

%% Eccentric Anomaly
n = sqrt(u./a.^3);

% eccentric anomaly
EA = unwrap(2*atan(((1-e)./(1+e)).^0.5.*tan(v/2)));

%% Time of Periapse Passage 
tp = t-1./n.*(EA-e.*sin(EA));















