% stateRates.m
% Author: Pierce Martin
% Calculates rates of change of state (including State Transition Matrix) 
% based on J2 and drag model
%
% Inputs:
%   t: current simulation time
%   state: Current simulation state in form [x y z xdot ydot zdot mu J2 CD
%   xs1 ys1 zs1 xs2 ys2 zs2 xs3 ys3 zs3 phi_1 .... phi_n]'
%
%   params: paramter file including RE, m, A, rho0, r0, H, theta_dot
%
% Outputs:
%   rates: rate of change of state vector
function [rates] = stateRates(t,state)

% Extract constants (for easier typing);
% RE = params.RE;
% m = params.m;
% A = params.A;
% rho0 = params.rho0;
% r0 = params.r0;
% H = params.H;
% theta_dot = params.theta_dot;

%% Define Constants
% Oblateness
% uE   = 3.986004415e14;                     % m^3/s^2
% J2   = 1.082626925638815e-3;                % []
RE  = 6378136.3;                           % m    Radius of earth

% Drag
% Cd      = 2.0;
A       = 3.0;                  % m^2
m       = 970;                  % kg
rho0    = 3.614e-13;            % kg/m^3
r0      = 700000.0 + RE;       % m
H       = 88667.0;              % m
theta_dot  = 7.29211585530066e-5;  % rad/s

% Extract state components
x = state(1,1);
y = state(2,1);
z = state(3,1);
xdot = state(4,1);
ydot = state(5,1);
zdot = state(6,1);
mu = state(7,1);
J2 = state(8,1);
CD = state(9,1);
xs1 = state(10,1);
ys1 = state(11,1);
zs1 = state(12,1);
xs2 = state(13,1);
ys2 = state(14,1);
zs2 = state(15,1);
xs3 = state(16,1);
ys3 = state(17,1);
zs3 = state(18,1);
phi = reshape(state(19:342,1),18,18);

% Calculate additional parameters for calculations
r     = sqrt(x^2+y^2+z^2);
rhoA  = rho0*exp(-(r-r0)/H);

VA    = sqrt((xdot+theta_dot*y)^2+(ydot-theta_dot*x)^2+zdot^2);
xddot = -mu*x/r^3*(1-3/2*J2*(RE/r)^2*(5*(z/r)^2-1))-1/2*CD*(A/m)*rhoA*VA*(xdot+theta_dot*y);
yddot = -mu*y/r^3*(1-3/2*J2*(RE/r)^2*(5*(z/r)^2-1))-1/2*CD*(A/m)*rhoA*VA*(ydot-theta_dot*x);
zddot = -mu*z/r^3*(1-3/2*J2*(RE/r)^2*(5*(z/r)^2-3))-1/2*CD*(A/m)*rhoA*VA*zdot;

% Calculate rate of change of X-state
F = transpose([xdot ydot zdot xddot yddot zddot 0 0 0 0 0 0 0 0 0 0 0 0]);

% Calculate A matirx
Amat = A_18x18(A,CD,H,J2,RE,m,r0,rho0,theta_dot,mu,x,xdot,y,ydot,z,zdot); 

% Calculate phi_dot
phi_dot = Amat*phi;

% Create rate vector
rates = [F;reshape(phi_dot,numel(phi_dot),1)];