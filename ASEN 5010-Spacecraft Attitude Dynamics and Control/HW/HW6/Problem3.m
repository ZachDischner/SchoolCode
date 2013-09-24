clc;clear all; close all

I           = diag( [10,20,30] );
sigma_init  = [0.8 0.1 -0.1]';
E321_init   = MRP2Euler321(sigma_init);
B_init      = MRP2EP(sigma_init);
w_init      = [0 2 0]';

L = [1 2 -1]';
%% Perform Integration for 321 Euler Angles
t = linspace(0,200,500);


%% Perform Integration for MRPs
K = 50;%1e-5*[1 0 0; 0 2 0; 0 0 3];  % Proportional gain matrix
P = 50;%1e-5*[1 0 0; 0 2 0; 0 0 3];  % Derivative Gain matrix
K1=.001;


z_init = dz(K,sigma_init);
X_init = [w_init; sigma_init; z_init]';

w(1,:) = w_init';
sigma(1,:) = sigma_init';
X(1,:) = X_init;




U_MRP = @(sigma,w,L,K,P,K1,z) -K*sigma - P*w - P*K*K1*z + P*K1*I*[w-[0 2 0]'];


for ii=2:length(t)
    %% Integrate
    % Find the Derivative
    Xprime = diffMRPw_Integral(t(ii-1),X(ii-1,:), I, L, U_MRP, K, P,K1);
    X(ii,:) = X(ii-1,:) + (t(ii) - t(ii-1))*Xprime';
    w(ii,:) = X(ii,1:3);
    sigma(ii,:) = X(ii,4:6);
end

% Plot the Angular Position
figure; set(gcf,'Color',[1 1 1], 'Position',[10 (900)  900 500])
subplot(2,1,1)
title('Regulator Problem With MRP Control Law');
hold on
plot(t,sigma(:,1),'r-'); plot(t,sigma(:,2),'k.'); plot(t,sigma(:,3),'b--');
xlabel('Time [s]'); ylabel('MRPs')
legend('$\sigma_1$','$\sigma_2$','$\sigma_3$','location','best')

% Plot the Angular Velocity
subplot(2,1,2)
hold on
plot(t,w(:,1),'r-'); plot(t,w(:,2),'k.'); plot(t,w(:,3),'b--');
xlabel('Time [s]'); ylabel('Angular Velocity [rad/s]')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')







