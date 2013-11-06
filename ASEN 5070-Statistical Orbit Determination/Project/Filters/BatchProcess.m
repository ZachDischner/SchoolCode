% function Xstar0 = BatchProcess()
% Compute new xhat
% Looks for functions to return:
%   G
%   Htilde
% Inside of /Filters/BatchTools.m

clc; clear all;close all; format compact;tic
warning off MATLAB:nearlySingularMatrix

%% First, Load and Parse Observation Data
load('Observations.mat');
t_obs       = obs(:,1);
station     = obs(:,2);
rho_obs     = obs(:,3);
rhodot_obs  = obs(:,4);

%% Pre-Allocations
x       = zeros(length(obs));
y       = x;
z       = x;
xdot    = x;
ydot    = x;
zdot    = x;
Xsite1  = x;
Ysite1  = x;
Zsite1  = x;
Xsite2  = x;
Ysite2  = x;
Zsite2  = x;
Xsite3  = x;
Ysite3  = x;
Zsite3  = x;
Phi     = cell(1,length(obs));
y_res   = zeros(2,length(obs));
Xsite   = x;
Ysite   = x;
Zsite   = x;
theta   = x;


%% Initialize Variables

% Calculations for rho, rhodot. Put into bigger functiom sometime
%---------------------------------------------
findrhostar      = @(x,y,z,Xsite,Ysite,Zsite,theta) sqrt(x^2+y^2+z^2+Xsite^2+Ysite^2+Zsite^2-2*(x*Xsite+y*Ysite)*cos(theta)+2*(x*Ysite-y*Xsite)*sin(theta)-2*z*Zsite);
findrhodotstar   = @(x,y,z,xdot,ydot,zdot,Xsite,Ysite,Zsite,theta,theta_dot,rho) (x*xdot + y*ydot + z*zdot - (xdot*Xsite + ydot*Ysite)*cos(theta) + theta_dot*(x*Xsite + y*Ysite)*sin(theta)...
                    +(xdot*Ysite - ydot*Xsite)*sin(theta) + theta_dot*(x*Ysite - y*Xsite)*cos(theta) - zdot*Zsite)...
                                                    /rho;
%---------------------------------------------

% System Constants
%---------------------------------------------
Phi_Init    = eye(18,18);
tol         = 1e-13;
uE          = 3.986004415e14;        % m^3/s^2
J2          = 1.082626925638815e-3;  % []
Cd          = 2;                     % []
theta_dot   = 7.29211585530066e-5;   % rad/s
time        = t_obs;
%---------------------------------------------

% Information
%---------------------------------------------
sigma_rho   = 0.01;                  % rms std
sigma_rhodot= 0.001;                 % rms std
R           = [sigma_rho^2  ,         0        ; ...
                    0       ,    sigma_rhodot^2];
W = inv(R); 
% Range only
% W(2,2)=0;
Pbar0       = diag([1e6,1e6,1e6,1e6,1e6,1e6,1e20,1e6,1e6,1e-10,1e-10,1e-10,1e6,1e6,1e6,1e6,1e6,1e6]);
%---------------------------------------------

% Initial Conditions
%---------------------------------------------
RV_Init     = [757700,5222607.0,4851500.0,2213.21,4678.34,-5371.30];
Station_Init= [-5127510.0 , -3794160.0 , 0.0 ,...               %101
                3860910.0 , 3238490.0  , 3898094.0 , ...        %337
                549505.0  , -1380872.0 , 6182197.0 ];           %394
Const_Init = [uE , J2 , Cd ];
%---------------------------------------------

% Form Initialization State
%---------------------------------------------
Xstar0 = [RV_Init , Const_Init , Station_Init , reshape(Phi_Init,1,length(Phi_Init)^2)]';
%---------------------------------------------

% Initial xbar0, or a-priori state deviation from reference trajectory
%---------------------------------------------
xbar0       = zeros(18,1);
%---------------------------------------------


%% Perform Batch Loop
num_iterations = 3;
res=[];
for ii = 1:num_iterations
    
    % Dynamical Integration
    %---------------------------------------------
    tol_mat     = ones(size(Xstar0)) .* tol;
    options     = odeset('RelTol',tol,'AbsTol',tol,'OutputFcn',@odetpbar);

    [time,StatePhi] = ode45(@StateDeriv_WithPhi,time,Xstar0,options);
    %---------------------------------------------
    
    
        % Batch Processing Part
    %---------------------------------------------
    Lam = inv(Pbar0);
    N   = Pbar0\xbar0;  % same as inv(pobar)*xbar0
    

    % Reform Phi Matrix
    %---------------------------------------------
    for jj = 1:length(time)
           Phi      = reshape(StatePhi(jj,19:end),size(Phi_Init));
           Xstar    = StatePhi(:,1:18);
           x        = Xstar(:,1);
           y        = Xstar(:,2);
           z        = Xstar(:,3);
           xdot     = Xstar(:,4);
           ydot     = Xstar(:,5);
           zdot     = Xstar(:,6);
           Xsite1   = Xstar(:,10);
           Ysite1   = Xstar(:,11);
           Zsite1   = Xstar(:,12);
           Xsite2   = Xstar(:,13);
           Ysite2   = Xstar(:,14);
           Zsite2   = Xstar(:,15);
           Xsite3   = Xstar(:,16);
           Ysite3   = Xstar(:,17);
           Zsite3   = Xstar(:,18);

    %---------------------------------------------
    
        theta(jj)          = theta_dot*time(jj);
        Htilde             = zeros(2,18);
        % Check Stations
        %---------------------------------------------
        %Station 1
        if station(jj) == 101
            Xsite(jj)   = Xsite1(jj);   Ysite(jj)=Ysite1(jj);   Zsite(jj)=Zsite1(jj);
            % Find H Tilde
            %---------------------------------------------
            Htilde = FindHtilde(Xsite(jj),Ysite(jj),Zsite(jj),theta(jj),theta_dot,x(jj),xdot(jj),y(jj),ydot(jj),z(jj),zdot(jj));
            %---------------------------------------------
            Htilde  = [Htilde , zeros(2,6)];
            
        end

        %Station 2
        if station(jj) == 337
            Xsite(jj)   = Xsite2(jj);   Ysite(jj)=Ysite2(jj);   Zsite(jj)=Zsite2(jj);
            % Find H Tilde
            %---------------------------------------------
            Htilde = FindHtilde(Xsite(jj),Ysite(jj),Zsite(jj),theta(jj),theta_dot,x(jj),xdot(jj),y(jj),ydot(jj),z(jj),zdot(jj));
            %---------------------------------------------
            Htilde  = [Htilde(:,1:9) , zeros(2,3), Htilde(:,10:12),zeros(2,3)];
        end

        %Station 3
        if station(jj) == 394
            Xsite(jj)   = Xsite3(jj);   Ysite(jj)=Ysite3(jj);   Zsite(jj)=Zsite3(jj);
            % Find H Tilde
            %---------------------------------------------
            Htilde = FindHtilde(Xsite(jj),Ysite(jj),Zsite(jj),theta(jj),theta_dot,x(jj),xdot(jj),y(jj),ydot(jj),z(jj),zdot(jj));
            %---------------------------------------------
            Htilde  = [Htilde(:,1:9),zeros(2,6),Htilde(:,10:12)];
        end
        %---------------------------------------------
       
        % Map To Epoch
        %---------------------------------------------
        H{jj}    = Htilde*Phi;
        %---------------------------------------------
        
        % Cumulate INFORMATION?? Matrix
        %---------------------------------------------
        Lam       = Lam + H{jj}'*W*H{jj};
        %---------------------------------------------
        
        % Put into FindG
        %---------------------------------------------
        rhostar   = findrhostar(x(jj),y(jj),z(jj),Xsite(jj),Ysite(jj),Zsite(jj),theta(jj));
        rhodotstar= findrhodotstar(x(jj),y(jj),z(jj),xdot(jj),ydot(jj),zdot(jj),Xsite(jj),Ysite(jj),Zsite(jj),theta(jj),theta_dot,rhostar);
        %---------------------------------------------
        
        % Find Observation Deviations
        %---------------------------------------------
        ystar       = [rhostar;rhodotstar];
        y_res(:,jj) = [rho_obs(jj);rhodot_obs(jj)] - ystar;
        %---------------------------------------------

        
        % Cumulate SOMETING?? Matrix
        %---------------------------------------------
        N     = N + H{jj}'*W*y_res(:,jj);
        %---------------------------------------------                                                         
    end
    fprintf('RMS of rho is :  %3.5f \n',rms(y_res(1,:)))
    fprintf('RMS of rhodot is :  %3.5f \n',rms(y_res(2,:)))
    % Find New State Deviation
    %---------------------------------------------
    xhat0   = Lam\(N);
    %---------------------------------------------
    
    % Update Best Guess of Initial Conditions
    %---------------------------------------------
    Xstar0  = [Xstar0(1:18) + xhat0; (reshape(Phi_Init,length(Phi_Init)^2,1))];
    %--------------------------------------------- 
    
    % Update a-priori State Deviation% function Xstar0 = BatchProcess()
    %---------------------------------------------
    xbar0 = xbar0 - xhat0;
    %---------------------------------------------
    
    
    figure(1)
    subplot(num_iterations,2,2*ii-1)
    plot(y_res(1,:))
    ylabel('$\rho$ residuals')
    xlabel('observation number')
    
    subplot(num_iterations,2,2*ii)
    plot(y_res(2,:))
    ylabel('$\dot{\rho}$ residuals')
    xlabel('observation number')
    
    res = [res,[rms(y_res(1,:));rms(y_res(2,:))]];
    
end

fprintf('\n\nRunning Time for Batch Processor : %3.5f\n\n',toc)

% For output stuff
% table=[Xstar0(1:18),diag(inv(Lam))];
Pfinal = Phi*inv(Lam)*Phi';
state=Xstar0(1:18);
cov     = diag(inv(Lam));
change = Xstar0(1:18)-[RV_Init , Const_Init , Station_Init]';
 save('Batch.mat','state','cov','Pfinal','change');

 % IF you want to save images
figure_awesome('save','eps')

% t=[Xstar0(1:18),Xstar0(1:18)-[RV_Init , Const_Init , Station_Init]', diag(inv(Lam))]












