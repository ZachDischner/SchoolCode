%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        getTropoCorrection.m
% Author      : Zach Dischner
% Date        : 10/30/2013
% Description : Retrieve tropospheic correction from simple model (in
%               meters)
% 
%
%                 __...____________________          ,
%                `(\ [ ===NCC-1700===--|__|) ___..--"_`--.._____
%                  `"""""""""""""""""| |""` [_""_-___________"_/
%                                    | |   /..../`'-._.-'`
%                                ____| |__/::..'_
%                               |\ ".`"` '_____//\
%                               `"'-.   """""  \\/
%                                    `""""""""""`
% Inputs        : ZenithCorr - Zenith delay [meters]
%                 el - elevation angle in degrees
% 
% Outputs       : TropoDelay - Delay value [meters]
% 
% History       October 30 2013 - First Rev
%               
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

function TropoDelay = getTropoCorrection(ZenithCorr, el)
TropoDelay = abs(ZenithCorr/sind(el));