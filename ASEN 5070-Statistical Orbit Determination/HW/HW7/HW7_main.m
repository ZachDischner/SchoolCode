%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Zach Dischner-10/21/2012
% 
% ASEN 5070-Statistical Orbit Determination
% 
% Homework 7
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;close all; format compact;format long g;tic



%% 1 - Derive A_18x18 and Htilde_2x18

% [A,Htilde]= FindA_Htilde();
% fprintf('A and Htilde were found. Not shown since they are too big to glean any insight from')

%% 2 - Perform integration
tol = 1e-12;
uE  = 398600.4 * 1000;          % m^3/s^2
J2  = 0.00108248;               % []
Cd  = 2;
theta_dot  = 7.29211585530066e-5;  % rad/s

time    = [0:20:18340];

RV_Init         = [757700,5222607.0,4851500.0,2213.21,4678.34,-5371.30];
Station_Init    = [-5127510.0 , -3794160.0 , 0.0 ,...               %101
                    3860910.0  , 3238490.0  , 3898094.0 , ...       %337
                    549505.0 , -1380872.0 , 6182197.0 ];            %394
Const_Init = [uE , J2 , Cd ];

Phi_Init = eye(18);



StateInit = [RV_Init , Const_Init , Station_Init , reshape(Phi_Init,1,length(Phi_Init)^2)]';

tol_mat = ones(size(StateInit)) .* tol;
options = odeset('RelTol',tol,'AbsTol',tol_mat);

[time,StatePhi] = ode45('StateDeriv_WithPhi',time,StateInit,options);

for ii = 1:length(time)
       Phi{ii} = reshape(StatePhi(ii,19:end),size(Phi_Init));
       State(ii,:) = StatePhi(ii,1:18);
end

fprintf('State(18340,0)  ==>   %f3.5 \n\n',State(time(end)/20,1))


%% 3 - Compare Residules of Observed and Computed Range and Range Rate 
rho = @(x,y,z,Xsite,Ysite,Zsite,theta) sqrt(x^2+y^2+z^2+Xsite^2+Ysite^2+Zsite^2-2*(x*Xsite+y*Ysite)*cos(theta)+2*(x*Ysite-y*Xsite)*sin(theta)-2*z*Zsite);
rhodot = @(x,y,z,xdot,ydot,zdot,Xsite,Ysite,Zsite,theta,theta_dot,rho) (x*xdot + y*ydot + z*zdot - (xdot*Xsite + ydot*Ysite)*cos(theta) + theta_dot*(x*Xsite + y*Ysite)*sin(theta)...
            +(xdot*Ysite - ydot*Xsite)*sin(theta) + theta_dot*(x*Ysite - y*Xsite)*cos(theta) - zdot*Zsite)...
                                                            /rho;
% Load in Observation Data
obs         = load('Observations.mat');
time_obs    = obs.obs(:,1);
station     = obs.obs(:,2);
rho_obs     = obs.obs(:,3);
rhodot_obs  = obs.obs(:,4);

SatState = cell(18,1);
rho_comp = zeros(length(time_obs),1); rhodot_comp = rho_comp;

for ii = 1:length(time_obs);
    
    index = time == time_obs(ii);
    SatState=num2cell(State(index,:));
    
    
    [x,y,z,xdot,ydot,zdot,uE,J2,Cd,Xsite1,Ysite1,Zsite1,Xsite2,Ysite2,Zsite2,Xsite3,Ysite3,Zsite3] = SatState{:};
    
    %Station 1
    if station(ii) == 101
        Xsite=Xsite1;   Ysite=Ysite1;   Zsite=Zsite1;
    end
    
    %Station 2
    if station(ii) == 337
        Xsite=Xsite2;   Ysite=Ysite2;   Zsite=Zsite2;
    end
    
    %Station 3
    if station(ii) == 394
        Xsite=Xsite3;   Ysite=Ysite3;   Zsite=Zsite3;
    end
    
    
    theta = (time_obs(ii)*theta_dot);
    rho_comp(ii)    = rho(x,y,z,Xsite,Ysite,Zsite,theta);
    rhodot_comp(ii) = rhodot(x,y,z,xdot,ydot,zdot,Xsite,Ysite,Zsite,theta,theta_dot,rho_comp(ii));
end

rho_diff = (rho_obs - rho_comp);
rhodot_diff = (rhodot_obs - rhodot_comp);

plot(time_obs,rho_diff)
ylabel('Rho Difference')
figure
plot(time_obs,rhodot_diff)
ylabel('Rho\dot Difference')

range_rms = sqrt(sum(rho_diff).^2/length(rho_diff));
range_rate_rms = sqrt(sum(rhodot_diff).^2/length(rhodot_diff));

fprintf('Range RMS is  ==>  %3.5f\n',range_rms)
fprintf('Range Rate RMS is  ==>  %3.5f\n\n',range_rate_rms)











fprintf('Time it took to run is : %f3.5',toc)
                