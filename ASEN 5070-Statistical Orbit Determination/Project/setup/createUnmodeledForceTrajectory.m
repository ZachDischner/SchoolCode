%% Setup
clc, close all, clear all

%% Load observation data
obs_data = load('ObservationData.txt');

% Extract data values
t_obs = obs_data(:,1);
station_id = obs_data(:,2);
rho_observed = obs_data(:,3);
rho_dot_observed = obs_data(:,4);

%% Known Parameters and values
params.RE = 6378136.3;
params.m = 970;
params.A = 3;
params.rho0 = 3.614e-13;
params.r0 = 700000.0+params.RE;
params.H = 88667.0;
params.theta_dot = 7.2921158553e-5;
ode_options = odeset('RelTol',1e-13,'AbsTol',1e-13,'OutputFcn',@odetpbar);
t0 = 0;

% Initial conditions
R0 = [757700.0;5222607.0;4851500.0];
V0 = [2213.21;4678.34;-5371.30];
mu0 = 3.986004415e14;
J2_0 = 1.082626925638815e-3;
CD_0 = 2.0;
XS1_0 = [-5127510.0;-3794160.0;0.0];
XS2_0 = [3860910.0;3238490.0;3898094.0];
XS3_0 = [549505.0;-1380872.0;6182197.0];
X_star0 = [R0;V0;mu0;J2_0;CD_0;XS1_0;XS2_0;XS3_0];

% Initialize data
rho_data = zeros(size(rho_observed));
rho_dot_data = rho_data;

%% Integrate trajectory with unmodeled force

% Create t_span to integrate over
t_data = [t_obs;t_obs+t_obs(end)+20];
t_data = [t_data;t_obs+t_data(end)+20];
t_data = [t_data;t_obs+t_data(end)+20];
station_id = repmat(station_id,4,1);

% Integrate state and phi
[tout state] = ode45(@stateRates_wUnmodeledForce,t_data,X_star0,ode_options,params);

% Exctract state components
x = state(:,1);
y = state(:,2);
z = state(:,3);
xdot = state(:,4);
ydot = state(:,5);
zdot = state(:,6);
mu = state(:,7);
J2 = state(:,8);
CD = state(:,9);
xs1 = state(:,10);
ys1 = state(:,11);
zs1 = state(:,12);
xs2 = state(:,13);
ys2 = state(:,14);
zs2 = state(:,15);
xs3 = state(:,16);
ys3 = state(:,17);
zs3 = state(:,18);

% Calculate range and range rate
for i = 1:length(t_data)
    % Calculate current theta value
    theta_dot = params.theta_dot;
    theta = theta_dot*(tout(i) - t0);
    
    % Determine which station the observation is from
    switch station_id(i)
        
        case 101
            xs = xs1(i);
            ys = ys1(i);
            zs = zs1(i);
            
        case 337
            xs = xs2(i);
            ys = ys2(i);
            zs = zs2(i);

            
        case 394
            xs = xs3(i);
            ys = ys3(i);
            zs = zs3(i);
    end
    
    % Create data with unmodeled force
    rho_data(i) = rho(theta,x(i),xs,y(i),ys,z(i),zs)+0*randn(1,1);
    rho_dot_data(i) = rho_dot(theta,theta_dot,x(i),xdot(i),xs,y(i),ydot(i),ys,z(i),zdot(i),zs)+0*randn(1,1);
end

% Create output data matrix
data_out = [t_data station_id rho_data rho_dot_data];
save('ObservationData_wUnmodeledForce.txt','data_out','-ascii', '-tabs');
