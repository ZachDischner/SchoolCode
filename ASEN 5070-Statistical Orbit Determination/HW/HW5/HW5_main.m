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
tol     = 1e-9;
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


%% 1c - Verify PHI Symplecticness from EQ 4.2.22

Phi_inv = [Phi{end}(3:end,3:end)',  -Phi{end}(1:2,3:end)';...
            -Phi{end}(3:end,1:2)',    Phi{end}(1:2,1:2)'];
fprintf('\n\n#^#^#^#^#^#^#^#^#^#^#^# Problem 1c #^#^#^#^#^#^#^#^#^#^#^#^#^#\n\n')
fprintf('Phi(t_100,t_0) * Phi(t_100,t_0)^-1 is nearly the identity matrix. Close enough for numerical precision\n')
Phi{end}*Phi_inv

% Could also do:
% J = [[0 0 ;0 0;-eye(2,2)],[eye(2,2);0 0;0 0]]
% Phi{end}*-(J*Phi{end}*J)'


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


















fprintf('\n\nAssignement took %3.3f seconds to run\n\n',toc)