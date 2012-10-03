%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Zach Dischner-8/28/2012
% 
% ASEN 5070-Statistical Orbit Determination
% 
% Homework 1
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath 'Utils'

clc;clear all;close all
tic

%% Problem 1-Convert Cartesian to Keplarian Coordinates
u   = 398600.5;     %km^3/s^2
R   = [-2436.45 -2436.45 6891.037];     % [i j k] km
V   = [5.088611 -5.088611 0.0];         % [i j k] km/s  => rDOT

[a e i Omega w nu EA T] = RVtoKepler(R,V,u)   % Gotta go in and look up that it is elliptical


%% Problem 2-Convert the Keplarian Coordinates Above Back to Cartesian
% [R V] = Keplar_to_RV(a,e,i,Omega,w,T); 
[R_new V_new] = Kepler_to_RV_easy(a,e,i,Omega,w,nu,u);

%% Problem 3-Solve for 2 body Acceleration Due to Gravity 
% A_x = -x*u/sqrt(x^2 + y^2 + z^2)^3    ...

%% Problem 4- Compute Future R,V of Spacecraft for 2 Orbit Periods
% Need 2 Full Orbit Periods
P       = 2*pi*sqrt((a^3)/u);

time    = [0:20:2*P]; 
tol     = 1e-19;
options = odeset('RelTol',tol,'AbsTol',[tol,tol,tol,tol,tol,tol]);

RV_init = [R';V'];

[t,RV]  = ode45('RV_Deriv1',time,RV_init,options);



%% Problem 5-Plot Radius, Velocity, and Accl vs Time 

figure

% Radius
r_norm = sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2);
subplot(3,1,1)
plot(t,r_norm,'linewidth',2)
title('Orbit Radius Vs Time '); xlabel('Time [s]'); ylabel('Radius [km]'); grid on

% Velocity
v_norm  = sqrt(RV(:,4).^2 + RV(:,5).^2 + RV(:,6).^2);
subplot(3,1,2)
plot(t,v_norm,'linewidth',2)
title('Orbit Velocity Vs Time '); xlabel('Time [s]'); ylabel('Velocity [km/s]'); grid on

% Acceleration
AV = RV_Deriv(-0,RV);
subplot(3,1,3)
plot(t,sqrt(AV(:,4).^2 + AV(:,5).^2 + AV(:,6).^2),'linewidth',2)
title('Orbit Acceleration Vs Time '); xlabel('Time [s]'); ylabel('Accleration [km/s^2]'); grid on


%% Problem 6-Compute Specific Kinetic and Potential NRG, Plot Change in Total Specific Energy 
% Kinetic
E_kinetic   = v_norm.^2/2;

% Potential 
E_potential = -u./r_norm;

% Total Specific Energy
TE          = E_kinetic + E_potential;

figure
plot(time,TE-TE(1))
xlabel('Time [s]'); ylabel('dTE'); title('Change in Total Energy Vs Time')


%% Problem 7-CH1 #1 From Text
X     = Problem7();
% fprintf(['X0 \t',num2str(X0),'\n'...
%          'Y0 \t',num2str(Y0),'\n'...
%          'Xp0 \t',num2str(Xp0),'\n'...
%          'Yp0 \t',num2str(Yp0),'\n'...
%          ]);





disp(['HW1 Took:   ',num2str(toc),'     seconds to run'])
% figure_awesome






