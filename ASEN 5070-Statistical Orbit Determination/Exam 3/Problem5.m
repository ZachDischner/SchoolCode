%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   ASEN 5070-Problem 5
% 
% Zach Dischner
%   Exam 3
%       Problem 5
% 
% 
% 
% Solves and answers questions relating to problem 5 of the STATOD final
% exam
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
clc;clear all;close all

sigma = 2.49;
Beta = 0.045;

%% Load Data
d=load('Exam3_Problem5_data.txt');
t = d(:,1);
ob = d(:,2);

Truth = sin(2*pi*t/10);

load('EtaHatEX3.mat')

%% Smooth out the Data
l = length(t);

% for kk = length(t)-1:-1:2
for ii = 1:998
   kk = l-ii;
   if ii == 1   
         
        S(kk) = P(kk)*PhiStep(kk)'*inv(P(kk+1));
        
        PS(kk) = P(kk) + S(kk)*P(kk+1)*S(kk)';
        
        EtaHatS(kk) = EtaHat(kk) + S(kk)*(EtaHat(kk)-PhiStep(kk)*EtaHat(kk));
   else
        S(kk) = P(kk)*PhiStep(kk)'*inv(P(kk+1));
        
        PS(kk) = P(kk) + S(kk)*(PS(kk+1) - P(kk+1))*S(kk)';
        
        EtaHatS(kk) = EtaHat(kk) + S(kk)*(EtaHat(kk+1) - PhiStep(kk+1)*EtaHat(kk));
   end   
end


figure
plot(Truth,'k','linewidth',1)
hold on
plot(ob,'r.','MarkerSize',1)
plot(EtaHat,'b','linewidth',1)
plot(EtaHatS,'g','Linewidth',1)
legend('Truth','Observations','$\hat{\eta}$','Smoothed $\hat{\eta}$')
xlabel('t');ylabel('$\hat{\eta}$');title('Comparison - Smoothing in Action')


%% Plot Residual Comparision
res = Truth-EtaHat';
res2 = Truth(1:999) - EtaHatS';
figure
hold on
plot(Truth,'k')
plot(Truth+res,'b','Linewidth',1)
plot(Truth(1:999)+res2,'g','Linewidth',1)

legend('Truth','Filtered RMS','Smoothed RMS')
xlabel('t');ylabel('$\eta$');title(' - Residual Comparision')



%% Normality Fit
figure
subplot(1,2,1)
hist(res);histfit(res)
title('Original Residuals')

subplot(1,2,2)
hist(res2);histfit(res2)
title('Smoothed Residuals')

%% RMS Calculation
fprintf('\n Original res: rms = %3.5f\n\n',rms(res))
fprintf('\n Smoothed res: rms = %3.5f\n\n',rms(res2))
fprintf('Resulting in a drop in res of:  %3.5f Percent\n\n',(rms(res)-rms(res2))*100)