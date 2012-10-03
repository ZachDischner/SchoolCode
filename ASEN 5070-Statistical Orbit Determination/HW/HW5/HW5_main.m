%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Zach Dischner-9/28/2012
% 
% ASEN 5070-Statistical Orbit Determination
% 
% Homework 5
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath 'Utils'

clc;clear all;close all; format compact;format short g;
% set(0,'defaultLineLineWidth', 2);set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on')
tic

%% Problem 1a-Simple Orbit Solver
%1a - Get true solution
RV_Init = [1.0 0 0 1.0]';
tol     = 1e-15;
tol_mat = ones(size(RV_Init)) .* tol;
time    = [0:0.01:100];

options = odeset('RelTol',tol,'AbsTol',tol_mat);


[t,X]  = ode45('RV_Deriv1a',time,RV_Init,options);

RV_true = X(mod(t,10)==0,:);

%% 1b Perturb previous initial conditions

perturb_init = 1e-6*[1 -1 1 1]';
RV_Init2    = RV_Init - perturb_init;

Phi_Init_mat= eye(size(RV_Init2,1));
[m n]       = size(Phi_Init_mat);

% Reform initial vector of phi. 
Phi_Init    = reshape(Phi_Init_mat,m*n,1);

State_Init  = [RV_Init2;Phi_Init];

tol_mat = ones(size(State_Init)) .* tol;
options = odeset('RelTol',tol,'AbsTol',tol_mat);

[t,Y]   = ode45('RV_Deriv1b',time,State_Init,options);

Xstar    = Y(:,1:length(RV_Init2));
R       = sqrt(Xstar(:,1).^2 + Xstar(:,2).^2);

for ii = 1:length(Y)
    Phi{ii} = reshape(Y(ii,5:20),size(Phi_Init_mat));
end


fprintf('Tolerance used was %3.3e\n',tol)
fprintf('\n\n#^#^#^#^#^#^#^#^#^#^#^# Problem 1b #^#^#^#^#^#^#^#^#^#^#^#^#^#\n\n')
disp(['X(t_100) :                    [',num2str(X(end,:)),']']) 
disp(['X*(t_100) :                   [',num2str(Xstar(end,:)),']'])
disp(['X*(t_100) - X*(t_100) :    [',num2str(X(end,:) - Xstar(end,:)),']'])
disp('\Phi(t_100) :')
Phi{end}


%% 1c - Verify PHI Symplecticness from EQ 4.2.22

Phi_inv = [Phi{end}(3:end,3:end)',  -Phi{end}(1:2,3:end)';...
            -Phi{end}(3:end,1:2)',    Phi{end}(1:2,1:2)'];
fprintf('\n\n#^#^#^#^#^#^#^#^#^#^#^# Problem 1c #^#^#^#^#^#^#^#^#^#^#^#^#^#\n\n')
fprintf('Phi(t_100,t_0) * Phi(t_100,t_0)^-1 is nearly the identity matrix. Close enough for numerical precision\n')
Phi{end}*Phi_inv


%% 1d - Calculate Perturbation Vector
% Method 1: DELTA(X)i = Xi - X*i

Perturb1 = Xstar - X;

% Method 2: DELTA(X)i = PHI(ti,t0)*(delta??)Xt0 
for ii = 1:length(Phi)
    Perturb2(ii,:) = Phi{ii}*perturb_init;
end

figure; grid on

subplot(3,1,1)
plot(time,Perturb1); xlabel('Time [s]');ylabel('\delta(X)_1');title('X* - X')
subplot(3,1,2)
plot(time,Perturb2); xlabel('Time [s]');ylabel('\delta(X)_2');title('\Phi * \delta_x')
subplot(3,1,3)
plot(time,Perturb2 - Perturb1); xlabel('Time [s]');ylabel('\delta(X)_2 - \delta(x)_1');title('Difference in Perturbation Methods')







%% Problem 2     given  y = H*x + e,   x is a scalar
% Observations
y       = [1 2 1]';
% Weighting
W       = [2 0 0;
           0 1 0;
           0 0 1];
% Map between observations and state
H       = [1 1 1]';
% Priori
x_bar   = 2;
W_bar   = 2;

%% 2a - Find "x" Using Batch Processor Algorithm?? Just single point>>

% Eq 4.3.25 in book. 
x_caret = inv(H'*W*H + W_bar)*(H'*W*y + W_bar*x_bar);

fprintf('\n\n\n\n#^#^#^#^#^#^#^#^#^#^#^# Problem 2a #^#^#^#^#^#^#^#^#^#^#^#^#^#\n\n')
fprintf('x^    : %3.5f\n\n',x_caret)
%% 2b - Find best estimate for observation error e^

% Eq 4.3.7
e_caret = y - H*x_caret;
fprintf('\n\n#^#^#^#^#^#^#^#^#^#^#^# Problem 2b #^#^#^#^#^#^#^#^#^#^#^#^#^#\n\n')
disp(['e^    :   [',num2str(e_caret'),']'''])



%% Problem 3 - Problem 4.15 From Book
rho = [     6.37687486186586
            5.50318198665912
            5.94513302809067
            6.30210798411686
            5.19084347133671
            6.31368240334678
            5.80399842220377
            5.45115048359871
            5.91089305965839
            5.67697312013520
            5.25263404969825];

rho_dot = [ -0.00317546143535849
             1.17587430814596
            -1.47058865193489
             0.489030779000695
             0.993054430595876
            -1.40470245576321
             0.939807575607138
             0.425908088320457
            -1.47604467619908
             1.42173765213734
            -0.12082311844776];
sigma_rho       = 0.25;
sigma_rhodot    = 0.1;

R       = [ 0.0625   0;    % R = E(e*e') = R = W^-1
               0   0.01;
          ];
      
      XV_Init = [3.0;0.0];       % [ m ; m/s ]   [x v]
      
      Phi_Init= eye(2,2);
      
          
          DelXbar0 = [0 ; 0];
      for jj = 1:4
         
          
          
          State_Init  = [XV_Init ; reshape(Phi_Init,4,1)];
          
          % Priori Reference Trajectory
          Xstar0  = [4.0 ; 0.2];
          
          % Priori <something>
          Pbar0   = [1000 0 ; 0 100];
          
          % System dynamics
          k1          = 2.5;  % N/m
          k2          = 3.7;  % N/m
          m           = 1.5;  % Kg
          w_squared   = (k1 + k2)/m;
          h           = 5.4;  % m
          
          
          %% Perform dynamical integration
          time    = [0:1:10];
          tol_mat = ones(size(State_Init)) .* tol;
          options = odeset('RelTol',tol,'AbsTol',tol_mat);
          [t,X]      = ode45('SpringDeriv',time,State_Init);
          
          Xstar = X(:,1:length(XV_Init));
          for ii = 1:length(X)
              Phi3{ii} = reshape(X(ii,3:6),size(Phi_Init));
          end
          
          x = X(:,1); v = X(:,2);
          sumHH = [0 0;0 0];
          sumHy = [0;0];
          
          
          for ii = 1:length(rho)
              Htilde   = [                 x(ii)/rho(ii)                   0;...
                  (v(ii)/rho(ii) - x(ii)^2*v(ii)/(rho(ii)^3))    x(ii)/rho(ii) ...
                  ];
              H3{ii}    = Htilde*Phi3{ii};
              
              sumHH    = sumHH + H3{ii}'*inv(R)*H3{ii};
              
              y3{ii}    = [rho(ii) rho_dot(ii)]' - [sqrt(x(ii)^2 + h^2);x(ii)*v(ii)/rho(ii)];
              sumHy    = sumHy + H3{ii}'*inv(R)*y3{ii};
              
          end
          
          DelXcarret0 = inv(sumHH + inv(Pbar0))*(sumHy + inv(Pbar0)*DelXbar0);
          
          DelXbar0 = DelXbar0 - DelXcarret0;
          
          XV_Init = XV_Init + DelXcarret0;
          
      end
      format long G
      
fprintf('\n\n\n\n#^#^#^#^#^#^#^#^#^#^#^# Problem 3 #^#^#^#^#^#^#^#^#^#^#^#^#^#\n\n')
disp(['X*0    :   [',num2str(XV_Init'),']'''])





















fprintf('\n\n\n\nAssignement took %3.3f seconds to run\n\n\n\n',toc)