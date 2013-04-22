%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           projectSims.m
% Author:   Zach Dischner
% Date:     April 6, 2013
% 
% Usage:
%   projectSims;
%
% Description:  Perform all project required simulations ad obtain all
% v             isuals requested
% 
% Inputs:  
%
% Outputs: 
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%% Prepare Workspace
clc;clear all;close all

%% Part A: Plot Satellite Position Vector
t   = linspace(0,10*60,500);
r_N = computeR_N(t);

figure;set(gcf,'Color',[1 1 1], 'Position',[10 (900)  1100 700])
plot(r_N);xlabel('Time [s]');ylabel('Orbit Position [km]');
legend('$^Nr_1$','$^Nr_2$','$^Nr_3$','location','best');
title('$^Nr$ Satellite Position Components as Functions of Time')


%% Part E: Simulate Magnetic Field Vector Evolution
M_N = zeros(size(r_N));
for ii=1:length(t)
    M_N(ii,:) = r2MagVec( r_N(ii,:), t(ii) ) ;
end

figure;set(gcf,'Color',[1 1 1], 'Position',[10 (900)  1100 700])
plot(M_N);xlabel('Time [s]');ylabel('Magnetic Field [nT]');
legend('$^NM_1$','$^NM_2$','$^NM_3$','location','best');
title('$^NM$ Earth''s Magnetic Field Components as Functions of Time')
