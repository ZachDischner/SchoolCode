%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Zach Dischner-9/8/2012
% 
% ASEN 5070-Statistical Orbit Determination
% 
% Homework 2
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath 'Utils'

clc;clear all;close all;format compact
tic

%% Problem 1-Integrate Equations of Motion With Oblateness

% 1A-Symbolic Derivative of Grad(U)

syms x y z r phi u R_E J2
U = (u/r)*(1-J2*(R_E/2)^2*(1.5*(z/r)^2-0.5));

U = subs(U,'r','sqrt(x^2 + y^2 + z^2)');


GradU = [diff(U,x); diff(U,y); diff(U,z)];
disp('Grad U is:')
pretty(GradU)

disp('A Simplified version is:')
GradU   = simplify(GradU);
pretty(subs(GradU,'sqrt(x^2 + y^2 + z^2)','r'))


% 1B-Plot Orbit Elements, Integrated for one day. 

u   = 398600.4;     % km^3/s^2
J2  = 0.00108248;   % []
R_E = 6378.1145;    % km    Radius of earth

% Initial conditions 
R   = [-2436.45 -2436.45 6891.0379];    % [i j k] km
V   = [5.088611 -5.088611 0.0];         % [i j k] km/s  => rDOT

time    = [0:20:24*60*60]; 
tol     = 1e-15;
options = odeset('RelTol',tol,'AbsTol',[tol,tol,tol,tol,tol,tol]);

RV_init = [R';V'];

% Perform Integration 
[t1,RV1]  = ode45('RV_Deriv_With_Oblateness',time,RV_init,options);
% save('RV_Computed.mat','t1','RV1');

% For fast debugging, just load saved version of t,RV
% load('RV_Computed.mat');

% Allocate and Convert Elements 
a1       = zeros(length(RV1),1);
e1       = a1;
i1       = a1;
Omega1   = a1;
w1       = a1;
nu1      = a1;
EA1      = a1;
T1       = a1;

[a1 e1 i1 Omega1 w1 nu1 Tp1] = MatRV2Kepler(RV1(:,1:3),RV1(:,4:6),u,t1);

%Do the plots 
figure
subplot(3,1,1)
plot(t1,a1,'linewidth',2);  xlabel('Time [s]'); ylabel('a [Km]');   title('a vs. time'); grid on

subplot(3,1,2)
plot(t1,e1,'linewidth',2);  xlabel('Time [s]'); ylabel('e [rad]');   title('Eccentricity vs. time'); grid on

subplot(3,1,3)
plot(t1,i1,'linewidth',2);  xlabel('Time [s]'); ylabel('i [rad]');   title('Inclination vs. time'); grid on

figure
subplot(3,1,1)
plot(t1,Omega1,'linewidth',2);  xlabel('Time [s]'); ylabel('\Omega [rad]');   title('RAAN vs. time'); grid on

subplot(3,1,2) % Fix this!!!
plot(t1,mod(w1,pi),'linewidth',2);  xlabel('Time [s]'); ylabel('w [rad]');   title('Argument of Periapse vs. time'); grid on

subplot(3,1,3) % Fix this!!!
plot(t1,Tp1,'linewidth',2);  xlabel('Time [s]'); ylabel('Tp');   title('Tp vs. time'); grid on


%   1C-Compute Specific Energy, verify constantness
x1 = RV1(:,1); y1 = RV1(:,2); z1 = RV1(:,3);
v1 = sqrt(RV1(:,4).^2 + RV1(:,5).^2 + RV1(:,6).^2);
r1 = sqrt(RV1(:,1).^2 + RV1(:,2).^2 + RV1(:,3).^2); 
U1 = (u./r1).*(1-J2*(R_E./r1).^2.*(3/2.*(z1./r1).^2-.5));
E1 = (v1.^2)./2 - U1;

figure
plot(t1,E1-E1(1),'linewidth',2)
xlabel('Time [s]'); ylabel('E - E_0 [km^2/s^2]'); title('Specific Energy vs Time With J2 Perturbation')



%   1D-Show Constant Angular Momentum
h1  = cross(RV1(:,1:3),RV1(:,4:6));
hk1 = h1(:,3);
figure
plot(t1,hk1 - hk1(1),'linewidth',2)
xlabel('Time [s]'); ylabel('h_k - h_k_0  [N*km*s]'); title('Angular Momentum h_k vs Time With J2 Perturbation')






%% Problem 2-Now Include Drag in Equations of Motion to Integrate 
[t2,RV2]  = ode45('RV_Deriv_With_Oblateness_Drag',time,RV_init,options);
save('RV_Computed_Drag.mat','t2','RV2');
% For fast debugging, just load saved version of t,RV
% load('RV_Computed_Drag.mat');


x2 = RV2(:,1); y2 = RV2(:,2); z2 = RV2(:,3);
v2 = sqrt(RV2(:,4).^2 + RV2(:,5).^2 + RV2(:,6).^2);
r2 = sqrt(RV2(:,1).^2 + RV2(:,2).^2 + RV2(:,3).^2); 
U2 = (u./r1).*(1-J2*(R_E./r1).^2.*(3/2.*(z1./r1).^2-.5));
E2 = (v2.^2)./2 - U2;

figure
plot(t2,E2-E2(1))
xlabel('Time [s]'); ylabel('E - E_0 [km^2/s^2]'); title('Specific Energy vs Time With J2 Perturbation and Drag')



% 2B-Changes in Orbital Elements from Drag
% Dumb way to do it
a2       = zeros(length(RV2),1);
e2       = a2;
i2       = a2;
Omega2   = a2;
w2       = a2;
nu2      = a2;
EA2      = a2;
T2       = a2;

[a2 e2 i2 Omega2 w2 nu2 Tp2] = MatRV2Kepler(RV2(:,1:3),RV2(:,4:6),u,t2);

figure
subplot(3,1,1)
plot(t1,a2-a1,'linewidth',2);  xlabel('Time [s]'); ylabel('\Delta a [Km]');   title('Change in a vs. time'); grid on

subplot(3,1,2)
plot(t1,e2-e1,'linewidth',2);  xlabel('Time [s]'); ylabel('\Delta e [rad]');   title('Change in Eccentricity vs. time'); grid on

subplot(3,1,3)
plot(t1,i2-i1,'linewidth',2);  xlabel('Time [s]'); ylabel('\Delta i [rad]');   title('Change in Inclination vs. time'); grid on

figure
subplot(3,1,1)
plot(t1,Omega2 - Omega1,'linewidth',2);  xlabel('Time [s]'); ylabel('\Delta \Omega [rad]');   title('Change in RAAN vs. time'); grid on

subplot(3,1,2) % Fix this!!!
plot(t1,w2 - w1,'linewidth',2);  xlabel('Time [s]'); ylabel('\Delta w [rad]');   title('Change in Argument of Periapse vs. time'); grid on

subplot(3,1,3) % Fix this!!!
plot(t1,Tp2 - Tp1,'linewidth',2);  xlabel('Time [s]'); ylabel('\Delta Tp');   title('Change in Tp vs. time'); grid on


% 



disp(['HW1 Took:   ',num2str(toc),'     seconds to run'])
figure_awesome('save');






