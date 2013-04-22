%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           intMRPW.m
% Author:   Zach Dischner
% Date:     April 20, 2013
% 
% Usage:
%   intMRPW;
%
% Description:  Numerically integrate attitude and omega. Project part f
% 
% Inputs:  
%
% Outputs: 
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
clc;clear all; close all;


% Variables 
%---------------------------------------------
% Initial Attitude
euler321_init    = deg2rad([5 10 -5])';     % Rad
% Initial Rate
w_init      = deg2rad([0.4 0.3 0.2])';      % Rad/s
I           = [ 25     2.5       0.5;...
                2.5     20        0;...
                0.5     0         15];
%---------------------------------------------

% Integration
%---------------------------------------------
t       = linspace(0,60*95,1000);%linspace(0,60*10,601);
X_init  = [w_init ; euler321_init];
X       = zeros(length(t),length(X_init));
X(1,:)  = X_init;
sigma   = zeros(length(t),length(euler321_init));
sigma(1,:) = MRPswitch( Euler3212MRP(euler321_init),1);
w       = zeros(length(t), length(w_init));
w(1,:)  = w_init;

% Initial Values at t=0
s_B     = zeros(length(t),3); 
s_B(1,:)= computeSunVec_B( sigma(1,:) );

M_B     = zeros(length(t),3);

r_N     = computeR_N(t(1));
EN      = ECI2ECF( t(1) );

r_ECEF  = EN*r_N';
M_N     = r2MagVec( r_ECEF, t(1) );


r_orbit = norm(r_N);
lam     = atan2(r_ECEF(2), r_ECEF(1)); 
phi     = asin(r_ECEF(3)/r_orbit);
                 

TE      = Earth2TopoDCM(lam,phi);

M_B(1,:)= TE*EN*M_N;
s_N     = [0 -1 0]'; % EN'*TE'*s_B(1,:)';
VB = [M_B(1,:)',s_B(1,:)'];
VN = [M_N,s_N];
sigma_OLAE(1,:)  = MRPswitch( Gibbs2MRP( OLAE(VB,VN)),1)';



for ii=2:length(t)
    
    %% Integrate
   % Find the Derivative
   Xprime = diff321w(t(ii-1),X(ii-1,:),I);
   
   % General linear integration:
   %    x_(n+1) = x_(n) + x'*delta_t
   X(ii,:) = X(ii-1,:) + Xprime'*(t(ii)-t(ii-1));
   sigma(ii,:) = MRPswitch( Euler3212MRP(X(ii,4:6)),1);
   w(ii,:) = X(ii,1:3);
   
   %% Get Sun and Magnetic Field Vectors
   s_B(ii,:)= computeSunVec_B( sigma(ii,:) );

   % Variable to hold Magnetic Field vector computations
   r_N     = computeR_N( t(ii-1) );
   r_orbit = norm( r_N );
   lam     = atan2(r_ECEF(2), r_ECEF(1)); % *look at it again
   phi     = asin(r_ECEF(3)/r_orbit);
   EN      = ECI2ECF( t(ii) );
   TE      = Earth2TopoDCM(lam,phi);
   
   r_ECEF  = EN*r_N';
   M_N     = r2MagVec( r_ECEF, t(ii) );
   M_B(ii,:)=TE*EN*M_N;
   s_N     = [0 -1 0]'; %EN'*TE'*s_B(ii,:)';   % N frame, s is always [0 -1 0]
   
   
%    sigmaTrue(ii,:)   = MRPswitch( C2MRP(TE*EN)' ,1 );
   VB = [M_B(ii,:)',s_B(ii,:)'];
   VN = [M_N,s_N];
   sigma_OLAE(ii,:)  = MRPswitch( Gibbs2MRP( OLAE(VB,VN)),1)';
%    sigmaDiff(:,ii) = sigmaTrue-sigma_OLAE;
   
%    sigma(ii,:) = MRPswitch(sigma(ii-1,:) + sigmaprime'*(t(ii)-t(ii-1)),1);

    %% Compare OLAE
end

%% MRP and Rate Simulation
figure;set(gcf,'Color',[1 1 1], 'Position',[10 (900)  900 500])
subplot(2,1,1)
plot(t,sigma)
title('MRP Simulation $\sigma_{B/N}$')
xlabel('Time [s]'); ylabel('Modified Rodreguizes')
legend('$\sigma_1$','$\sigma_2$','$\sigma_3$','location','best')

subplot(2,1,2)
plot(t,w);title('Rate Simulation $w_{B/N}$')
xlabel('Time [s]'); ylabel('[Rad/S]')
legend('$w_1$','$w_2$','$w_3$','location','best')

%% Field Vector Sets
figure;set(gcf,'Color',[1 1 1], 'Position',[10 (900)  900 500])
subplot(2,1,1)
plot(t,s_B)
title('Simulated Sun Vector $^Bs$')
xlabel('Time [s]'); ylabel('Sun Vector Components')
legend('$^Bs_1$','$^Bs_2$','$^Bs_3$','location','best')

subplot(2,1,2)
plot(t,M_B)
title('Simulated Magnetic $^BM$')
xlabel('Time [s]'); ylabel('Magnetic Field Vector Components')
legend('$^BM_1$','$^BM_2$','$^BM_3$','location','best')


%% Difference in Attitudes
figure;set(gcf,'Color',[1 1 1], 'Position',[10 (900)  900 500])
subplot(3,1,1)
plot(t,sigma_OLAE(:,1),t,sigma(:,1))
title('Attitude Estimation Algorithms')
xlabel('Time [s]'); ylabel('$\sigma_1$')
legend('OLAE','True','location','best')

subplot(3,1,2)
plot(t,sigma_OLAE(:,2),t,sigma(:,2))
title('Attitude Estimation Algorithms')
xlabel('Time [s]'); ylabel('$\sigma_2$')
legend('OLAE','True','location','best')

subplot(3,1,3)
plot(t,sigma_OLAE(:,3),t,sigma(:,3))
title('Attitude Estimation Algorithms')
xlabel('Time [s]'); ylabel('$\sigma_3$')
legend('OLAE','True','location','best')






