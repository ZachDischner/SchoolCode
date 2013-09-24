clc;clear all; close all

I           = diag( [10,20,30] );
sigma_init  = [0.8 0.1 -0.1]';

w_init      = [0 0 0]';

L = [0 0 0]';
%% Perform Integration for 321 Euler Angles
t = linspace(0,1000,100);


% %% Perform Integration for MRPs
% K = 50;%1e-5*[1 0 0; 0 2 0; 0 0 3];  % Proportional gain matrix
% P = 50;%1e-5*[1 0 0; 0 2 0; 0 0 3];  % Derivative Gain matrix
% K1=.001;


% z_init = dz(K,sigma_init);
X_init = [sigma_init; w_init]';

w(1,:) = w_init';
sigma(1,:) = sigma_init';
X(1,:) = X_init;

P = diag([0.2, 0.28284, 0.3464]);
K = 0.004;


% U_MRP = @(sigma,w,L,K,P,K1,z) -K*sigma - P*w - P*K*K1*z + P*K1*I*[w-[0 2 0]'];


for ii=2:length(t)
    %% Integrate
    % Find the Derivative
    Xprime = diffMRPLin(t,X(ii-1,:),I,P,K);
    X(ii,:) = X(ii-1,:) + (t(ii) - t(ii-1))*Xprime';
    w(ii,:) = X(ii,4:6);
    sigma(ii,:) = X(ii,1:3);
end

% Plot the Angular Position
figure; set(gcf,'Color',[1 1 1], 'Position',[10 (900)  900 500])
subplot(2,2,1)
title('Linearized CLD');
hold on
plot(t,sigma(:,1),'r-'); plot(t,sigma(:,2),'k.'); plot(t,sigma(:,3),'b--');
xlabel('Time [s]'); ylabel('MRPs')
legend('$\sigma_1$','$\sigma_2$','$\sigma_3$','location','best')

% Plot the Angular Velocity
subplot(2,2,3)
hold on
plot(t,w(:,1),'r-'); plot(t,w(:,2),'k.'); plot(t,w(:,3),'b--');
xlabel('Time [s]'); ylabel('Angular Velocity [rad/s]')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')



%% Nonlinear CLD
clear X w sigma Xprime ii
sigma(1,:) = sigma_init;
w(1,:) = w_init;
X(1,:) = [sigma_init;w_init]';
% X_init = [w_init; sigma_init]';
% 
% w(1,:) = w_init';
% sigma(1,:) = sigma_init';
% X(1,:) = X_init;
% 
% 
% U_MRP = @(sigma,w,L,K,P) -K*sigma - P*w;
% 
% 
% for ii=2:length(t)
%     %% Integrate
%     % Find the Derivative
%     Xprime = diffMRPw(t(ii-1),X(ii-1,:), I, L, U_MRP, K, P);
%     X(ii,:) = X(ii-1,:) + (t(ii) - t(ii-1))*Xprime';
%     w(ii,:) = X(ii,1:3);
%     sigma(ii,:) = X(ii,4:6);
% end

for jj=2:length(t)
    %% Integrate
    % Find the Derivative
    Xprime = diffMRPw_part4(t,X(jj-1,:),I,P,K);
    X(jj,:) = X(jj-1,:) + (t(jj) - t(jj-1))*Xprime';
    w(jj,:) = X(jj,4:6);
    sigma(jj,:) = X(jj,1:3);
end

subplot(2,2,2)
title('Nonlinear CLD');
hold on
plot(t,sigma(:,1),'r-'); plot(t,sigma(:,2),'k.'); plot(t,sigma(:,3),'b--');
xlabel('Time [s]'); ylabel('MRPs')
legend('$\sigma_1$','$\sigma_2$','$\sigma_3$','location','best')

% Plot the Angular Velocity
subplot(2,2,4)
hold on
plot(t,w(:,1),'r-'); plot(t,w(:,2),'k.'); plot(t,w(:,3),'b--');
xlabel('Time [s]'); ylabel('Angular Velocity [rad/s]')
legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')




