%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Zach Dischner-12/4/2012
% 
% ASEN 5070-Statistical Orbit Determination
% 
% Homework 11
% 
% Inputs    : None
% 
% Outputs   : Plots for homework 11
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;close all; format compact;format long g;tic

%% Plot the Ellipsoid
tmp = load('P.mat');
P_Batch = tmp.P;
P = P_Batch(1:3,1:3);
[evecs,evals] = eig(P);
semi(1) = sqrt(evals(1,1));
semi(2) = sqrt(evals(2,2));
semi(3) = sqrt(evals(3,3));
semi = sort(semi);
semi = semi([3,2,1]);
plotEllipsoid(evecs,semi)


%% Problem 4-41, call function
Problem4_41

figure_awesome('save')