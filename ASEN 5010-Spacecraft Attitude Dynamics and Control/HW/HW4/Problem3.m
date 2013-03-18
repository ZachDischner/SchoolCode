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
t = linspace(0,100,1000);
I = [125  0   0;
     0   100  0;
     0    0  75];

w1_init = [1;0;0];
w1=intRigidEOM(w1_init,t,I);

w2_init = [0;1;0];
w2=intRigidEOM(w2_init,t,I);

w3_init = [0;0;1];
w3=intRigidEOM(w3_init,t,I);



figure;set(gcf,'Color',[1 1 1], 'Position',[10 (900)  900 500])
subplot(3,2,1)
plot(w1);xlabel('Time [s]');ylabel('$\omega$ [rad/s]');
title('Pure Spin About $\omega_1$')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')

subplot(3,2,3)
plot(w2);xlabel('Time [s]');ylabel('$\omega$ [rad/s]');
title('Pure Spin About $\omega_2$')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')

subplot(3,2,5)
plot(w3);xlabel('Time [s]');ylabel('$\omega$ [rad/s]');
title('Pure Spin About $\omega_3$')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')



%% Now with perturbations in other axes
w11_init = [1;0.1;0];
w11=intRigidEOM(w11_init,t,I);

w22_init = [0.1;1;0];
w22=intRigidEOM(w22_init,t,I);

w33_init = [0;0.1;1];
w33=intRigidEOM(w33_init,t,I);


subplot(3,2,2)
plot(w11);xlabel('Time [s]');ylabel('$\omega$ [rad/s]');
title('Spin About $\omega_1$ with $\omega_2$ Perturb')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')

subplot(3,2,4)
plot(w22);xlabel('Time [s]');ylabel('$\omega$ [rad/s]');
title('Spin About $\omega_2$ with $\omega_1$ Perturb')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')

subplot(3,2,6)
plot(w33);xlabel('Time [s]');ylabel('$\omega$ [rad/s]');
title('Spin About $\omega_3$ with $\omega_2$ Perturb')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')





