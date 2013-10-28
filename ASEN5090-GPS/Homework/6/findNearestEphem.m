%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        findNearestEmph.m
% Author      : Zach Dischner
% Date        : 10/11/2013
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
% Inputs        : PRN - Satellite PRN number
%                 GPSWeeks - GPS week number (modded or no?)
%                 GPSSOW - GPS Seconds of week
%                 navData - A full array of all emphimeris data, fetched
%                            from navigation file
% Outputs       : emphData - Single row (struct?) of emphemeris data per
%                            sat PRN at time [gps_weeks, gps_seconds
% 
% History       Oct 11 2013 - First Version
%               Oct 22 2013 - Added return for rownums
%               Oct 24 2013 - Changed PRN matching to ismember(), to allow
%                               for array matching of PRNs
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

function [ephemData,rownums] = findNearestEphem(PRN, GPS_Weeks, GPS_SOW, navData)

% weeknums = nav_ephem(:,19);
% secofweeks = nav_ephem(:,17)

rownums = find(navData(:,17)<=GPS_SOW & ismember(navData(:,1),PRN) & navData(:,19)==GPS_Weeks);
ephemData = navData(rownums,:);

