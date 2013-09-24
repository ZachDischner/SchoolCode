clc;clear all; close all

I           = diag( [10,20,30] );
sigma_init  = [0 0.8 0]';
E321_init   = MRP2Euler321(sigma_init);
B_init      = MRP2EP(sigma_init);
w_init      = [0 2 0]';

for kk = 1:2
    
    if kk == 1
        L = [0 0 0]';
    else
        L = [1 2 -1]';
    end
    %% Perform Integration for 321 Euler Angles
    t = linspace(0,15);
    X_init = [w_init; E321_init]';
    
    w(1,:) = w_init';
    ang(1,:) = E321_init';
    clear X
    X(1,:) = X_init;
    
    
    K = 50;%1e-5*[1 0 0; 0 2 0; 0 0 3];  % Proportional gain matrix
    P = 50;%1e-5*[1 0 0; 0 2 0; 0 0 3];  % Derivative Gain matrix
    
    % Good check of pos def
    % positivedefinite = all(eig(A) > 0);
    
    U_321 = @(theta,w,L,K,P) -(theta'*K*BmatEuler321(theta))' - P*w;
    
    
    for ii=2:length(t)
        %% Integrate
        % Find the Derivative
        Xprime = diff321w(t(ii-1), X(ii-1,:), I, L, U_321, K, P);
        X(ii,:) = X(ii-1,:) + (t(ii) - t(ii-1))*Xprime';
        w(ii,:) = X(ii,1:3);
        ang(ii,:) = X(ii,4:6);
    end
    
    % Plot the Angular Position
    figure; set(gcf,'Color',[1 1 1], 'Position',[10 (900)  900 500])
    subplot(2,1,1)
    title('Regulator Problem With Euler Angle Control Law');
    hold on
    plot(t,ang(:,1),'r-'); plot(t,ang(:,2),'k.'); plot(t,ang(:,3),'b--');
    xlabel('Time [s]'); ylabel('Euler Angles (321)')
    legend('$\theta_1$','$\theta_2$','$\theta_3$','location','best')
    
    % Plot the Angular Velocity
    subplot(2,1,2)
    hold on
    plot(t,w(:,1),'r-'); plot(t,w(:,2),'k.'); plot(t,w(:,3),'b--');
    xlabel('Time [s]'); ylabel('Angular Velocity [rad/s]')
    legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')
    
    if kk == 2
        fprintf('For the Euler angle 321 case, with unmodeled torques of:\n')
        fprintf('L = \t [%3.5f %3.5f %3.5f]\n\n',L(1),L(2),L(3))
        predicted = L'/K*BmatEuler321([0 0 0]);
        fprintf('Predicted theta_ss \t\t= [%3.5f %3.5f %3.5f]',predicted(1),predicted(2),predicted(3))
        act=ang(end,:);
        fprintf('\nAcutal simulated theta_ss \t= [%3.5f %3.5f %3.5f]\n\n',act(1),act(2),act(3))
    end
    
    
    %% Perform Integration for MRPs
    X_init = [w_init; sigma_init]';
    
    w(1,:) = w_init';
    sigma(1,:) = sigma_init';
    X(1,:) = X_init;
    
    
    U_MRP = @(sigma,w,L,K,P) -K*sigma - P*w;
    
    
    for ii=2:length(t)
        %% Integrate
        % Find the Derivative
        Xprime = diffMRPw(t(ii-1),X(ii-1,:), I, L, U_MRP, K, P);
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
    
        
    if kk == 2
        fprintf('\n\nFor the MRP set case, with unmodeled torques of:\n')
        fprintf('L = \t [%3.5f %3.5f %3.5f]\n\n',L(1),L(2),L(3))
        predicted = L'/K;
        fprintf('Predicted theta_ss \t\t= [%3.5f %3.5f %3.5f]',predicted(1),predicted(2),predicted(3))
        act=sigma(end,:);
        fprintf('\nAcutal simulated theta_ss \t= [%3.5f %3.5f %3.5f]\n\n',act(1),act(2),act(3))
    end
    
    
    
    %% Perform Integration for Euler Parameters
    X_init = [w_init; B_init]';
    
    w(1,:) = w_init';
    B(1,:) = B_init';
    clear X;
    X(1,:) = X_init;
    
    
    U_EP = @(epsilon,w,L,K,P) -K*epsilon - P*w;
    
    
    for ii=2:length(t)
        %% Integrate
        % Find the Derivative
        Xprime = diffQw(t(ii-1),X(ii-1,:), I, L, U_EP, K, P);
        X(ii,:) = X(ii-1,:) + (t(ii) - t(ii-1))*Xprime';
        w(ii,:) = X(ii,1:3);
        B(ii,:) = X(ii,4:7);
    end
    
    % Plot the Angular Position
    figure; set(gcf,'Color',[1 1 1], 'Position',[10 (900)  900 500])
    subplot(2,1,1)
    title('Regulator Problem With Euler Parameter Control Law');
    hold on
    plot(t,B(:,2),'r-'); plot(t,B(:,3),'k.'); plot(t,B(:,4),'b--');
    xlabel('Time [s]'); ylabel('MRPs')
    legend('$\epsilon_1$','$\epsilon_2$','$\epsilon_3$','location','best')
    
    % Plot the Angular Velocity
    subplot(2,1,2)
    hold on
    plot(t,w(:,1),'r-'); plot(t,w(:,2),'k.'); plot(t,w(:,3),'b--');
    xlabel('Time [s]'); ylabel('Angular Velocity [rad/s]')
    legend('$\omega_1$','$\omega_2$','$\omega_3$','location','best')
    
    if kk == 2
        fprintf('\n\nFor the Euler angle 321 case, with unmodeled torques of:\n')
        fprintf('L = \t [%3.5f %3.5f %3.5f]\n\n',L(1),L(2),L(3))
        predicted = L'/K;
        fprintf('Predicted theta_ss \t\t= [%3.5f %3.5f %3.5f]',predicted(1),predicted(2),predicted(3))
        act=B(end,2:4);
        fprintf('\nAcutal simulated theta_ss \t= [%3.5f %3.5f %3.5f]\n\n',act(1),act(2),act(3))
    end
    
end




