%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        findFirstEpoch.m
% Author      : Zach Dischner
% Date        : 10/24/2013
% Description : Function to return all emphimeris data from a nav data
%               array
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
% Inputs        : navData - Navigation dataset. 
% Outputs       : emphData - rows (struct?) of emphemeris data for
%                            the first epoch
%                 rows - row indices of the first epoch datasets
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

function [emphData,rows] = findFirstEpoch( navData )

weeknums = navData(:,19);
secofweeks = navData(:,17);

n_epochs = length(navData);
epochs = zeros(n_epochs,1);
for ii =1:n_epochs
    epochs(ii) = datenum(GPSTime2Date(weeknums(ii),secofweeks(ii)));
end

rows    = find(epochs==min(epochs));
emphData = navData(rows,:);