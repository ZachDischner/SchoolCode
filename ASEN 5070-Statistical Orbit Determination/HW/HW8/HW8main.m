%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Zach Dischner-10/31/2012
% 
% ASEN 5070-Statistical Orbit Determination
% 
% Homework 8
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;close all; format compact;format long g;tic



%% 1 - Different values of P

% Look in my function!!!

%% 2 - Plot for different values of P
e           = logspace(-15,-6,1000);
P2_True     = zeros(size(e));
P2_Kalman   = P2_True;
P2_Joseph   = P2_True;
P2_Potter   = P2_True;
P2_Batch    = P2_True;

for ii = 1:length(e)
    [P2_True(ii),P2_Kalman(ii),P2_Joseph(ii),P2_Potter(ii),P2_Batch(ii)]=FindP2(e(ii),1);
end

figure
subplot(4,1,1)
loglog(e,abs(P2_True-P2_Kalman));
xlabel('\epsilon');ylabel('P2_{true} - P2_{Kalman}'); title('Kalman')
subplot(4,1,2)
loglog(e,abs(P2_True-P2_Joseph));
xlabel('\epsilon');ylabel('P2_{true} - P2_{Joseph}'); title('Joseph')
subplot(4,1,3)
loglog(e,abs(P2_True-P2_Potter));
xlabel('\epsilon');ylabel('P2_{true} - P2_{Potter}'); title('Potter')
subplot(4,1,4)
loglog(e,abs(P2_True-P2_Batch));
xlabel('\epsilon');ylabel('P2_{true} - P2_{Batch}'); title('Batch')                   





%% 4
x_batch = zeros(size(e));
x_kalman = x_batch;
x_joseph = x_batch;
x_potter = x_batch;

for ii = 1:length(e)
    [x_batch(ii),x_kalman(ii),x_joseph(ii),x_potter(ii)] = FindStateP4(e,ii);
end
syms e
R       = eye(2,2);
H       = [1 2*e;1 3*e];
P0bar   = [1/e^2 0;0 1/e^2];
A       = inv(P0bar) + transpose(H)*R*H;
Ptrue   = inv(A);

















