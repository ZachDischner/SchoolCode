%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Zach Dischner-10/14/2012
% 
% ASEN 5070-Statistical Orbit Determination
% 
% Homework 6
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;close all; format compact;format long g;


%% Problem 4.20

[a b c] = Problem420();
fprintf('4.20-a           xhat  =  %3.5f \n\n',a)
fprintf('4.20-b    std(error)   =  %3.5f \n\n',b)
disp(['4.20-c        error    = [',num2str(c'),']'''])
fprintf('\n\n')


%% Problem 4.21
format rat
[a b c] = Problem421();
format rat
disp(['4.21-a        Batch Xhat0    = [',num2str(a'),']'''])
disp(['4.21-b        Sequential Xhat0    = [',num2str(b'),']'''])
disp(['4.21-c        Phi(1,0)*Xhat1 -->Xhat0    = [',num2str(c'),']'''])

fprintf('\n\n Or in a more RATIONAL way...\n')
a
b
c