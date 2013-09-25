%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        HW3.m
% Author      : Zach Dischner
% Date        : 9/24/2013
% Description : Matlat script for all calculations required for
%               ASEN 5090 Homework 3
%
%                 __...____________________          ,
%                `(\ [ ===NCC-1700===--|__|) ___..--"_`--.._____
%                  `"""""""""""""""""| |""` [_""_-___________"_/
%                                    | |   /..../`'-._.-'`
%                                ____| |__/::..'_
%                               |\ ".`"` '_____//\
%                               `"'-.   """""  \\/
%                                    `""""""""""`
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Setup Work Space
clc;clear all;close all
screen_size = get(0,'ScreenSize');
sw = screen_size(3);    % Screen Width
sh = screen_size(4);    % Screen Height
figColor = [0.99 0.99 0.98];
addpath gps_vis_ASEN5090_2013/