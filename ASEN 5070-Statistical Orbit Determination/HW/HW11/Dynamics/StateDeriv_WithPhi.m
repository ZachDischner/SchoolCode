function Xprime = StateDeriv_WithPhi(time,X)
%
%Author: 
%           Zach Dischner
%Date:
%           10/18/2012
%Use: 
%           StateHistory=ode45(StateDeriv_WithPhi,[linspace(1,100),X_init)
% 
% Derivative function, designed to return the derivative of the Cartesian
%   R and V observations, incliuding the effect of the Earth's Oblateness
%   and drag. In addition, return the PHI matrix, reformed into a single
%   vector. 
%
% Written by Zach Dischner on 8/30/2012
%
% inputs:
%         time: time...
%         State  : A column vector. 18-state + PHI reformed
%                       X
%                       Y
%                       Z
%                       X'
%                       Y'
%                       Z'
%                       u_e
%                       J2
%                       Cd
%                       Xbar_sitei   (1 2 or 3) based on station coordinates
%                       Phi
%
%

%% Define Constants
R_E     = 6378136.3;                           % m    Radius of earth
Area    = 3.0;                  % m^2
m       = 970;                  % kg
rho0    = 3.614e-13;            % kg/m^3
r0      = 700000.0 + R_E;       % m
H       = 88667.0;              % m
theta_dot  = 7.29211585530066e-5;  % rad/s


%% Form State Variables
Phi=reshape(X(19:end),sqrt(length(X(19:end))),sqrt(length(X(19:end))));
x       = X(1,1);
y       = X(2,1);
z       = X(3,1);
xdot    = X(4,1);
ydot    = X(5,1);
zdot    = X(6,1);
uE      = X(7,1);
J2      = X(8,1);
Cd      = X(9,1);
% Xsite1  = X(10,1);
% Ysite1  = X(11,1);
% Zsite1  = X(12,1);
% Xsite2  = X(13,1);
% Ysite2  = X(14,1);
% Zsite2  = X(15,1);
% Xsite3  = X(16,1);
% Ysite3  = X(17,1);
% Zsite3  = X(18,1);
r=sqrt(x^2 + y^2 + z^2);


%% Use Force Equations
F_U     = [...
              -uE/r^3*x*(1-3/2*J2*(R_E/r)^2*(5*(z/r)^2-1));
              -uE/r^3*y*(1-3/2*J2*(R_E/r)^2*(5*(z/r)^2-1));
              -uE/r^3*z*(1-3/2*J2*(R_E/r)^2*(5*(z/r)^2-3));
          ];
      
% [xdot ydot zdot] due to atmospheric drag
Va      = [
    xdot + theta_dot*y;
    ydot - theta_dot*x;
    zdot
    ];
va      = sqrt((xdot + theta_dot*y)^2 + (ydot - theta_dot*x)^2 + zdot^2);

rho_a   = rho0.*exp(-(r - r0)./H);

F_Drag = -0.5 .* Cd .* (Area./m) .* rho_a .* va .* Va;

% Assemble
F_a = F_U + F_Drag;

State_Deriv = [xdot ; ydot ; zdot ; F_a ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0]; % Add in station location with earth's rotation rate?

%% Integrate Phi
 A = FindA(Area,Cd,H,J2,R_E,m,r0,rho0,theta_dot,uE,x,xdot,y,ydot,z,zdot); 
bounds = length(Phi);
PhiPrime = A(1:bounds,1:bounds)*Phi;

%% Reform full state vector
Xprime = [State_Deriv ; reshape(PhiPrime,length(PhiPrime)^2,1)];











end











