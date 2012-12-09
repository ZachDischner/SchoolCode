      % function Xstar0 = ExtendedKalmanFilter()
% Compute new xhat
% Looks for functions to return:
%   G
%   Htilde
% Inside of /Filters/KalmanTools.m
%
% Notation:    t0 is a priori, like M_(t_i-1)
%              t1 is current step   M_(ti)

clc; clear all;close all; format compact;tic
warning off MATLAB:nearlySingularMatrix
set(0,'defaulttextinterpreter','latex')

prekalman = 1;



%% First, Load and Parse Observation Data
load('Observations.mat');
t_obs       = obs(:,1);
station     = obs(:,2); station(end+1) = station(end);station(end+1) = station(end);station(end+1) = station(end);
rho_obs     = obs(:,3);
rhodot_obs  = obs(:,4);

%% Pre-Allocations
x       = zeros(length(obs),1);
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
y1   = zeros(2,length(obs));
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
tol         = 1e-9;
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



if prekalman == 1
    num_iterations=3;
    cutoff_time=180*60;
    output=0;
    Xstar0=KalmanFunction(num_iterations,cutoff_time,output);
end
%     %% Start Preprocessing Kalman Loop
%     cutoff_time    = 2*180*60; % minutes
%     
%     t=time(time<cutoff_time);t(t==0)=[];t=[0;t];
%     
%     % Initialize Kalman Filter
%     %---------------------------------------------
%     xbar       = zeros(18,1);
%     xhat       = zeros(18,length(t));
%     %---------------------------------------------
%     
%     
%     %% Perform Kalman Loop
%     num_iterations = 3;
%     for ii = 1:num_iterations
%         tr=[];
%         P(:,:,1)          = Pbar0;
%         
%         % Dynamical Integration
%         %---------------------------------------------
%         tol_mat     = ones(size(Xstar0)) .* tol;
%         options     = odeset('RelTol',tol,'AbsTol',tol,'OutputFcn',@odetpbar);
%         
%         [t,StatePhi] = ode45(@StateDeriv_WithPhi,t,Xstar0,options);
%         %---------------------------------------------
%         
%         
%         % Reform Phi Matrix
%         %---------------------------------------------
%         for jj = 1:length(t)
%             Phi(:,:,jj)     = reshape(StatePhi(jj,19:end),size(Phi_Init));
%             Xstar    = StatePhi(:,1:18);
%             x        = Xstar(:,1);
%             y        = Xstar(:,2);
%             z        = Xstar(:,3);
%             xdot     = Xstar(:,4);
%             ydot     = Xstar(:,5);
%             zdot     = Xstar(:,6);
%             Xsite1   = Xstar(:,10);
%             Ysite1   = Xstar(:,11);
%             Zsite1   = Xstar(:,12);
%             Xsite2   = Xstar(:,13);
%             Ysite2   = Xstar(:,14);
%             Zsite2   = Xstar(:,15);
%             Xsite3   = Xstar(:,16);
%             Ysite3   = Xstar(:,17);
%             Zsite3   = Xstar(:,18);
%             
%             %---------------------------------------------
%             
%             
%             % Htilde, yi, Ki
%             %---------------------------------------------
%             theta(jj)          = theta_dot*time(jj);
%             Htilde             = zeros(2,18);
%             % Check Stations
%             %---------------------------------------------
%             %Station 1
%             if station(jj) == 101
%                 Xsite(jj)   = Xsite1(jj);   Ysite(jj)=Ysite1(jj);   Zsite(jj)=Zsite1(jj);
%                 % Find H Tilde
%                 %---------------------------------------------
%                 Htilde = FindHtilde(Xsite(jj),Ysite(jj),Zsite(jj),theta(jj),theta_dot,x(jj),xdot(jj),y(jj),ydot(jj),z(jj),zdot(jj));
%                 %---------------------------------------------
%                 Htilde  = [Htilde , zeros(2,6)];
%                 
%             end
%             
%             %Station 2
%             if station(jj) == 337
%                 Xsite(jj)   = Xsite2(jj);   Ysite(jj)=Ysite2(jj);   Zsite(jj)=Zsite2(jj);
%                 % Find H Tilde
%                 %---------------------------------------------
%                 Htilde = FindHtilde(Xsite(jj),Ysite(jj),Zsite(jj),theta(jj),theta_dot,x(jj),xdot(jj),y(jj),ydot(jj),z(jj),zdot(jj));
%                 %---------------------------------------------
%                 Htilde  = [Htilde(:,1:9) , zeros(2,3), Htilde(:,10:12),zeros(2,3)];
%             end
%             
%             %Station 3
%             if station(jj) == 394
%                 Xsite(jj)   = Xsite3(jj);   Ysite(jj)=Ysite3(jj);   Zsite(jj)=Zsite3(jj);
%                 % Find H Tilde
%                 %---------------------------------------------
%                 Htilde = FindHtilde(Xsite(jj),Ysite(jj),Zsite(jj),theta(jj),theta_dot,x(jj),xdot(jj),y(jj),ydot(jj),z(jj),zdot(jj));
%                 %---------------------------------------------
%                 Htilde  = [Htilde(:,1:9),zeros(2,6),Htilde(:,10:12)];
%             end
%             %---------------------------------------------
%             
%             
%             
%             % Time Update
%             %---------------------------------------------
%             if jj > 1
%                 Phi_step= Phi(:,:,jj)/Phi(:,:,jj-1);
%                 xbar(:,jj)   = Phi_step*xhat(:,jj-1);
%                 P(:,:,jj)   = Phi_step*P(:,:,jj-1)*Phi_step';
%                 
%             end
%             %---------------------------------------------
%             
%             % Put into FindG
%             %---------------------------------------------
%             rhostar   = findrhostar(x(jj),y(jj),z(jj),Xsite(jj),Ysite(jj),Zsite(jj),theta(jj));
%             rhodotstar= findrhodotstar(x(jj),y(jj),z(jj),xdot(jj),ydot(jj),zdot(jj),Xsite(jj),Ysite(jj),Zsite(jj),theta(jj),theta_dot,rhostar);
%             %---------------------------------------------
%             
%             % Find Observation Deviations
%             %---------------------------------------------
%             ystar       = [rhostar;rhodotstar];
%             y1(:,jj)   = [rho_obs(jj);rhodot_obs(jj)] - ystar;
%             %---------------------------------------------
%             
%             
%             % Kalman Gain
%             %---------------------------------------------
%             K1          = P(:,:,jj)*Htilde'*inv(Htilde*P(:,:,jj)*Htilde' + R);
%             %---------------------------------------------
%             
%             % Measurement Update
%             %---------------------------------------------
%             xhat(:,jj) = xbar(:,jj) + K1*(y1(:,jj) - Htilde*xbar(:,jj));
%             %         P(:,:,jj) = (eye(size(K1*Htilde)) - K1*Htilde)*P(:,:,jj)*(eye(size(K1*Htilde))-K1*Htilde)' + K1*R*K1';
%             P(:,:,jj) = (eye(size(K1*Htilde)) - K1*Htilde)*P(:,:,jj);
%             %---------------------------------------------
%             
%             tr = [tr,trace(P(1:3,1:3,jj))];
%             
%         end   % End observation loop
%         
%         
%         
%         % Update Best Guess of Initial Conditions
%         %---------------------------------------------%
%         Xstar0 = [Xstar0(1:18) + inv(Phi(:,:,end))*xhat(:,end); (reshape(Phi_Init,length(Phi_Init)^2,1))];
%         P      = inv(Phi(:,:,end))*P(:,:,end)*inv(Phi(:,:,end))';
%         xbar   = xbar(:,end-1) - xhat(:,end-1);
%         
%     end
% end


 cutoff_time=5*60;
 % Initialize Kalman Filter
%---------------------------------------------
xbar       = zeros(18,1);
xhat       = zeros(18,length(time));
%---------------------------------------------
num_iterations = 1;

textprogressbar('EKF Progress : ');
 for ii = 1:num_iterations
    tr=[];
    P(:,:,1)          = Pbar0;
    
    % Dynamical Integration Tolerances
    %---------------------------------------------
    tol_mat     = ones(size(Xstar0)) .* tol;
    options     = odeset('RelTol',tol,'AbsTol',tol);
    %---------------------------------------------
    
    
    StatePhi = Xstar0;

    StSw = [0,0];
    for jj = 1:length(time)
        progress=100 * jj/length(time);
        textprogressbar(progress);
        
        % Perform Dynamical Integration
        %---------------------------------------------
        if jj > 1
            [t,tmp] = ode45(@StateDeriv_WithPhi,[time(jj-1),time(jj)],StatePhi(:,jj-1),options);
            StatePhi(:,jj) = tmp(end,:)';
        else
            StatePhi(:,jj) = Xstar0;
        end
        %---------------------------------------------
           
        % Reform Phi Matrix
        %---------------------------------------------   
        Phi(:,:,jj)  = reshape(StatePhi(19:end,jj),size(Phi_Init));
        Xstar        = StatePhi(1:18,jj);
        x(jj)        = Xstar(1);
        y(jj)        = Xstar(2);
        z(jj)        = Xstar(3);
        xdot(jj)     = Xstar(4);
        ydot(jj)     = Xstar(5);
        zdot(jj)     = Xstar(6);
        Xsite1(jj)   = Xstar(10);
        Ysite1(jj)   = Xstar(11);
        Zsite1(jj)   = Xstar(12);
        Xsite2(jj)   = Xstar(13);
        Ysite2(jj)   = Xstar(14);
        Zsite2(jj)   = Xstar(15);
        Xsite3(jj)   = Xstar(16);
        Ysite3(jj)   = Xstar(17);
        Zsite3(jj)   = Xstar(18);
        
        %---------------------------------------------
    
    
        % Htilde, yi, Ki
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
        
        
        
        % Time Update
        %---------------------------------------------
        if jj > 1
%              Phi_step= Phi(:,:,jj)/Phi(:,:,jj-1);
            Phi_step    = Phi(:,:,jj);
            xbar(:,jj)  = Phi_step*xhat(:,jj-1);
            P(:,:,jj)   = Phi_step*P(:,:,jj-1)*Phi_step';
        end

        
        % Check Station Switchover
        if jj > 2
            if station(jj) ~= station(jj+1) || station(jj) ~= station(jj-1) || station(jj) ~= station(jj+2)
                StSw = [StSw;jj,1];
                factor = .001;
                P(:,:,jj)=P(:,:,jj).*factor;
            end
        end
            
        
        
        %---------------------------------------------
            
        % Put into FindG
        %---------------------------------------------
        rhostar   = findrhostar(x(jj),y(jj),z(jj),Xsite(jj),Ysite(jj),Zsite(jj),theta(jj));
        rhodotstar= findrhodotstar(x(jj),y(jj),z(jj),xdot(jj),ydot(jj),zdot(jj),Xsite(jj),Ysite(jj),Zsite(jj),theta(jj),theta_dot,rhostar);
        %---------------------------------------------
        
        % Find Observation Deviations
        %---------------------------------------------
        ystar       = [rhostar;rhodotstar];
        y1(:,jj)   = [rho_obs(jj);rhodot_obs(jj)] - ystar;
        %---------------------------------------------
        
        
        % Kalman Gain
        %---------------------------------------------
        K1          = P(:,:,jj)*Htilde'*inv(Htilde*P(:,:,jj)*Htilde' + R);
        %---------------------------------------------
        
        % Measurement Update
        %---------------------------------------------
        if time(jj) < cutoff_time
            xhat(:,jj) = zeros(18,1);% xhat(:,:,jj) =  xhat(:,:,jj); %xbar(:,:,jj) + K1*(y1(:,jj) - Htilde*xbar(:,:,jj));
%             xhat(:,jj) = xbar(:,jj) + K1*(y1(:,jj) - Htilde*xbar(:,jj));
        else
%             if ii < num_iterations
%                 % Update Best Guess of Initial Conditions
%                 %---------------------------------------------
%                 Xstar0  = [Xstar0(1:18) + inv(Phi_map)*xhat(:,jj-1);(reshape(Phi_map,length(Phi_Init)^2,1))];
%                 P      = inv(Phi_map)*P(:,:,end)*inv(Phi_map)';
%                 xbar   = xbar(:,end-1) - xhat(:,jj-1);
%                 %---------------------------------------------
%                 break
%             end
            if max(StSw(:,1) == jj)==1
                P(:,:,jj)=P(:,:,jj)./factor;                
            end
            xhat(:,jj) = K1*y1(:,jj);
         end
%         xhat(:,:,jj) = xbar(:,:,jj) + K1*(y1(:,jj) - Htilde*xbar(:,:,jj));
%         P(:,:,jj) = (eye(size(K1*Htilde)) - K1*Htilde)*P(:,:,jj)*(eye(size(K1*Htilde))-K1*Htilde)' + K1*R*K1';
        P(:,:,jj) = (eye(size(K1*Htilde)) - K1*Htilde)*P(:,:,jj);
        %---------------------------------------------
        
        % Update Nominal Trajectory
        %---------------------------------------------
        % Update StatePhi? Or make a State PhiBar??
        StatePhi(:,jj)  = [Xstar + xhat(:,jj); reshape(eye(size(Phi(:,:,jj))),length(Phi_Init)^2,1)];%(reshape(Phi(:,:,jj),length(Phi_Init)^2,1))];
        %---------------------------------------------
        

        tr = [tr,trace(P(1:3,1:3,jj))];
 
                                                               
    end   % End observation loop
    
    % Get rid of 0,0 element in station switch array
    StSw(1,:) = [];
        
    
    fprintf('\n\nRMS of rho is :  %3.5f \n',rms(y1(1,:)))
    fprintf('RMS of rhodot is :  %3.5f \n',rms(y1(2,:)))
    rho_comp = y1(1,:); rhodot_comp=y1(2,:);
    rho_comp(rho_comp<mean(rho_comp)-2*std(rho_comp) | rho_comp>mean(rho_comp)+2*std(rho_comp)) = [];
    rhodot_comp(rhodot_comp<mean(rhodot_comp)-2*std(rhodot_comp) | rhodot_comp>mean(rhodot_comp)+2*std(rhodot_comp)) = [];
    fprintf('\n\nWithout Outliers RMS of rho is :  %3.5f \n',rms(rho_comp))
    fprintf('Without Outliers RMS of rhodot is :  %3.5f \n',rms(rhodot_comp))

    
    y1_sw = y1(1,StSw(:,1))';
    y2_sw = y1(2,StSw(:,1))';
    
    for kk=1:length(StSw(:,1))
        if StSw(kk,1)> length(rho_comp)
            StSw(kk,:) = [0,0];
        end
    end
    StSw(StSw==0)=[];
    
    rho_sw = rho_comp(1,StSw(:,1))';
    rhodot_sw = rhodot_comp(1,StSw(:,1))';
    
    figure(1)
    subplot(num_iterations,2,2*ii-1)
    plot(y1(1,:))
    hold on
    plot(StSw(:,1),y1_sw,'r.','MarkerSize',20)
    ylabel('$\rho$ residules')
    xlabel('Observation Number')
    legend('Residules','Station Switchover')
    
    subplot(num_iterations,2,2*ii)
    plot(y1(2,:))
    hold on
    plot(StSw(:,1),y2_sw,'r.','MarkerSize',20)
    ylabel('$\dot{\rho}$ residules')
    xlabel('Observation Number')
    legend('Residules','Station Switchover')
    
    
    figure(2)
    subplot(num_iterations,1,ii)
    semilogy(tr)
    xlabel('Observation Number')
    ylabel('Trace( $P_{xyz}$ )')
    
    figure(3)
    subplot(num_iterations,2,2*ii-1)
    plot(rho_comp(1,:))
    hold on
    plot(StSw(:,1),rho_sw,'r.','MarkerSize',20)
    ylabel('$\rho$ residules')
    xlabel('Observation Number')
    legend('Residules','Station Switchover')
    
    subplot(num_iterations,2,2*ii)
    plot(rhodot_comp)
    hold on
    plot(StSw(:,1),rhodot_sw,'r.','MarkerSize',20)
    ylabel('$\dot{\rho}$ residules')
    xlabel('Observation Number')
    legend('Residules','Station Switchover')
    
 end
 
 
 figure
 plot(time(time<cutoff_time)/60,rho_obs(time<cutoff_time),'*r')
 hold on
 plot(time(time>cutoff_time)/60,rho_obs(time>cutoff_time),'.b')
 xlabel('time');ylabel('$\rho_{obs}$')
 legend('CKF','EKF Switch')

fprintf('\n\nRunning Time for Kalman Filter : %3.5f\n\n',toc)


















