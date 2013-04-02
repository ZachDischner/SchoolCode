%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           IntRigidEOM.m 
% Author:   Zach Dischner
% Date:     March 19, 2013
% 
% Usage:
%   
%
% Description:  
% 
% Inputs:  
%
% Outputs: 
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Setup
clc;clear all;close all
t = linspace(0,50,10000);
I = [125  0   0;
     0   100  0;
     0    0  75];
 
ang_init = [0;0;0];

w1_init = [1;0;0];
X1_init = [w1_init;ang_init];

X1=intRigidEOM(X1_init,t,I);

w2_init = [0;1;0];
X2_init = [w2_init;ang_init];
X2=intRigidEOM(X2_init,t,I);

w3_init = [0;0;1];
X3_init = [w3_init;ang_init];
X3=intRigidEOM(X3_init,t,I);



figure;set(gcf,'Color',[1 1 1], 'Position',[10 (900)  1100 700])
title('321 Euler Angle Sequence Attitudes')
subplot(3,2,1)
plotZ(t,X1(:,4:6));xlabel('Time [s]');ylabel('Attitude (321) [rad]');
title('Pure Spin About $\omega_1$')
legend('$\psi$','$\phi$','$\theta$','location','best')

subplot(3,2,3)
plotZ(t,X2(:,4:6));xlabel('Time [s]');ylabel('Attitude (321) [rad]');
title('Pure Spin About $\omega_2$')
legend('$\psi$','$\phi$','$\theta$','location','best')

subplot(3,2,5)
plotZ(t,X3(:,4:6));xlabel('Time [s]');ylabel('Attitude (321) [rad]');
title('Pure Spin About $\omega_3$')
legend('$\psi$','$\phi$','$\theta$','location','best')



%% Now with perturbations in other axes
w11_init = [1;0.1;0];
X11_init = [w11_init;ang_init];
X11=intRigidEOM(X11_init,t,I);

w22_init = [0.1;1;0];
X22_init = [w22_init;ang_init];
X22=intRigidEOM(X22_init,t,I);

w33_init = [0;0.1;1];
X33_init = [w33_init;ang_init];
X33=intRigidEOM(X33_init,t,I);


subplot(3,2,2)
plotZ(t,X11(:,4:6));xlabel('Time [s]');ylabel('Attitude (321) [rad]');
title('Spin About $\omega_1$ with $\omega_2$ Perturb')
legend('$\psi$','$\phi$','$\theta$','location','best')

subplot(3,2,4)
plotZ(t,X22(:,4:6));xlabel('Time [s]');ylabel('Attitude (321) [rad]');
title('Spin About $\omega_2$ with $\omega_1$ Perturb')
legend('$\psi$','$\phi$','$\theta$','location','best')

subplot(3,2,6)
plotZ(t,X33(:,4:6));xlabel('Time [s]');ylabel('Attitude (321) [rad]');
title('Spin About $\omega_3$ with $\omega_2$ Perturb')
legend('$\psi$','$\phi$','$\theta$','location','best')




figure;set(gcf,'Color',[1 1 1], 'Position',[10 (900)  1100 700])
title('Angular Velocity Components')
subplot(3,2,1)
plotZ(t,X1(:,1:3));xlabel('Time [s]');ylabel('Angular Velocies [rad/s]');
title('Pure Spin About $\omega_1$')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')

subplot(3,2,3)
plotZ(t,X2(:,1:3));xlabel('Time [s]');ylabel('Angular Velocies [rad/s]');
title('Pure Spin About $\omega_2$')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')

subplot(3,2,5)
plotZ(t,X3(:,1:3));xlabel('Time [s]');ylabel('Angular Velocies [rad/s]');
title('Pure Spin About $\omega_3$')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')

subplot(3,2,2)
plotZ(t,X11(:,1:3));xlabel('Time [s]');ylabel('Angular Velocies [rad/s]');
title('Spin About $\omega_1$ with $\omega_2$ Perturb')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')

subplot(3,2,4)
plotZ(t,X22(:,1:3));xlabel('Time [s]');ylabel('Angular Velocies [rad/s]');
title('Spin About $\omega_2$ with $\omega_1$ Perturb')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')

subplot(3,2,6)
plotZ(t,X33(:,1:3));xlabel('Time [s]');ylabel('Angular Velocies [rad/s]');
title('Spin About $\omega_3$ with $\omega_2$ Perturb')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')



%% Problem 4-Duffing 
X11=X22;
% Calculate duffing constants
I1=I(1,1); I2 = I(2,2); I3 = I(3,3);
w1 = X11(:,1); w2 = X11(:,2); w3 = X11(:,3);
H = sqrt(((I1.*w1).^2 + (I2.*w2).^2 + (I3.*w3).^2)); % for this case, w2 = w3 = 0


T = 1/2*I1*w1.^2 + 1/2*I2*w2.^2 + 1/2*I2*w2.^2;

A1 = (( I1-I2 ).*( 2*I3.*T-H.^2 ) + ( I1-I3 ).*( 2.*I2.*T-H.^2))./(I1*I2*I3);
A2 = (( I2-I3 ).*( 2*I1.*T-H.^2 ) + ( I2-I1 ).*( 2.*I3.*T-H.^2))./(I1*I2*I3);
A3 = (( I3-I1 ).*( 2*I2.*T-H.^2 ) + ( I3-I2 ).*( 2.*I1.*T-H.^2))./(I1*I2*I3);

B1 = 2*( I1-I2 )*( I1-I3 )/(I2*I3)*ones(length(A1),1);
B2 = 2*( I2-I1 )*( I2-I3 )/(I1*I3)*ones(length(A1),1);
B3 = 2*( I3-I2 )*( I3-I2 )/(I1*I2)*ones(length(A1),1);

returnPrime=1;
Xprime = intRigidEOM(X11_init,t,I,returnPrime);
w1dot=Xprime(:,1);
w2dot=Xprime(:,2);
w3dot=Xprime(:,3);

K1 = w1dot.^2 + A1.*w1.^2 + B1./2.*w1.^4;
K2 = w2dot.^2 + A2.*w2.^2 + B2./2.*w2.^4;
K3 = w3dot.^2 + A3.*w3.^2 + B3./2.*w3.^4;

K1true = (( 2*I2.*T - H.^2 ).*( H.^2 - 2.*I3.*T ))./(I1.^2.*I2.*I3);
K2true = (( 2*I3.*T - H.^2 ).*( H.^2 - 2.*I1.*T ))./(I2.^2.*I1.*I3);
K3true = (( 2*I1.*T - H.^2 ).*( H.^2 - 2.*I2.*T ))./(I3.^2.*I2.*I3);

fprintf('\n\nFor spin about b1 with perturbations about b2, the max difference in K values is:\n')
fprintf('K1 : [%3.5f]\n',max(K1-K1true))
fprintf('K2 : [%3.5f]\n',max(K2-K3true))
fprintf('K3 : [%3.5f]\n',max(K3-K3true))

figure
subplot(3,1,1);plot(t,K1-K1true);xlabel('Time [s]');ylabel('$dK_1$');
title('Verification of "energy type" Integral For Perturbed b1 Spin')
subplot(3,1,2);plot(t,K2-K2true);xlabel('Time [s]');ylabel('$dK_2$');
subplot(3,1,3);plot(t,K3-K3true);xlabel('Time [s]');ylabel('$dK_3$');

% pure spin about 2
% Calculate duffing constants
I1=I(1,1); I2 = I(2,2); I3 = I(3,3);
w1 = X1(:,1); w2 = X1(:,2); w3 = X1(:,3);
H = sqrt(((I1.*w1).^2 + (I2.*w2).^2 + (I3.*w3).^2)); % for this case, w2 = w3 = 0


T = 1/2*I1*w1.^2 + 1/2*I2*w2.^2 + 1/2*I2*w2.^2;

A1 = (( I1-I2 ).*( 2*I3.*T-H.^2 ) + ( I1-I3 ).*( 2.*I2.*T-H.^2))./(I1*I2*I3);
A2 = (( I2-I3 ).*( 2*I1.*T-H.^2 ) + ( I2-I1 ).*( 2.*I3.*T-H.^2))./(I1*I2*I3);
A3 = (( I3-I1 ).*( 2*I2.*T-H.^2 ) + ( I3-I2 ).*( 2.*I1.*T-H.^2))./(I1*I2*I3);

B1 = 2*( I1-I2 )*( I1-I3 )/(I2*I3)*ones(length(A1),1);
B2 = 2*( I2-I1 )*( I2-I3 )/(I1*I3)*ones(length(A1),1);
B3 = 2*( I3-I2 )*( I3-I2 )/(I1*I2)*ones(length(A1),1);

returnPrime=1;
Xprime = intRigidEOM(X1_init,t,I,returnPrime);
w1dot=Xprime(:,1);
w2dot=Xprime(:,2);
w3dot=Xprime(:,3);

K1 = w1dot.^2 + A1.*w1.^2 + B1./2.*w1.^4;
K2 = w2dot.^2 + A2.*w2.^2 + B2./2.*w2.^4;
K3 = w3dot.^2 + A3.*w3.^2 + B3./2.*w3.^4;

K1true = (( 2*I2.*T - H.^2 ).*( H.^2 - 2.*I3.*T ))./(I1.^2.*I2.*I3);
K2true = (( 2*I3.*T - H.^2 ).*( H.^2 - 2.*I1.*T ))./(I2.^2.*I1.*I3);
K3true = (( 2*I1.*T - H.^2 ).*( H.^2 - 2.*I2.*T ))./(I3.^2.*I2.*I3);
fprintf('\n\nFor pure spin about b1, the max difference in K values is:\n')
fprintf('K1 : [%3.5f]\n',max(K1-K1true))
fprintf('K2 : [%3.5f]\n',max(K2-K3true))
fprintf('K3 : [%3.5f]\n',max(K3-K3true))

figure
subplot(3,1,1);plot(t,K1-K1true);xlabel('Time [s]');ylabel('$dK_1$');
title('Verification of "energy type" Integral for Pure b1 spin')
subplot(3,1,2);plot(t,K2-K2true);xlabel('Time [s]');ylabel('$dK_2$');
subplot(3,1,3);plot(t,K3-K3true);xlabel('Time [s]');ylabel('$dK_3$');

% % Pure spin about 3
% I1=I(1,1); I2 = I(2,2); I3 = I(3,3);
% w1 = X1(:,1); w2 = X1(:,2); w3 = X1(:,3);
% H = I1.*w3; % for this case, w2 = w3 = 0
% T = 1/2*I1*w3.^2;
% 
% A1 = (( I1-I2 ).*( 2*I3.*T-H.^2 ) + ( I1-I3 ).*( 2.*I2.*T-H.^2))./(I1*I2*I3);
% A2 = (( I2-I3 ).*( 2*I1.*T-H.^2 ) + ( I2-I1 ).*( 2.*I3.*T-H.^2))./(I1*I2*I3);
% A3 = (( I3-I1 ).*( 2*I2.*T-H.^2 ) + ( I3-I2 ).*( 2.*I1.*T-H.^2))./(I1*I2*I3);
% 
% B1 = 2*( I1-I2 )*( I1-I3 )/(I2*I3)*ones(length(A1),1);
% B2 = 2*( I2-I1 )*( I2-I3 )/(I1*I3)*ones(length(A1),1);
% B3 = 2*( I3-I2 )*( I3-I2 )/(I1*I2)*ones(length(A1),1);
% 
% returnPrime=1;
% X1prime = intRigidEOM(X3_init,t,I,returnPrime);
% w1dot=X1prime(:,1);
% w2dot=X1prime(:,2);
% w3dot=X1prime(:,3);
% 
% K1 = w1dot.^2 + A1.*w1.^2 + B1./2.*w1.^4;
% K2 = w2dot.^2 + A2.*w2.^2 + B2./2.*w2.^4;
% K3 = w3dot.^2 + A3.*w3.^2 + B3./2.*w3.^4;
% 
% K1true = (( 2*I2.*T - H.^2 ).*( H.^2 - 2.*I3.*T ))./(I1.^2.*I2.*I3);
% K2true = (( 2*I3.*T - H.^2 ).*( H.^2 - 2.*I1.*T ))./(I2.^2.*I1.*I3);
% K3true = (( 2*I1.*T - H.^2 ).*( H.^2 - 2.*I2.*T ))./(I3.^2.*I2.*I3);
% 
% fprintf('\n\nFor pure spin about b3, the max difference in K values is:\n')
% fprintf('K1 : [%3.5f]\n',max(K1-K1true))
% fprintf('K2 : [%3.5f]\n',max(K2-K3true))
% fprintf('K3 : [%3.5f]\n',max(K3-K3true))




