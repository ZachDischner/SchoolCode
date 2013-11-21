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
%               Oct 31 2013 - Changed all terms to datenums, to account for
%                               week changeover
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

function [ephemData,rownums] = findNearestEphem(PRN, GPS_Weeks, GPS_SOW, navData)

% weeknums = nav_ephem(:,19);
% secofweeks = nav_ephem(:,17)
% sec_diff = abs(navData(:,17)-GPS_SOW);
% rownums = find( (sec_diff) == min(sec_diff) & ismember(navData(:,1),PRN) & navData(:,19)==GPS_Weeks);
% rownums = find( navData(:,17)<=GPS_SOW & ismember(navData(:,1),PRN) & navData(:,19)==GPS_Weeks);
GPSNUMBOOL = 1;
epoch_time  = GPSTime2Date(GPS_Weeks, GPS_SOW, GPSNUMBOOL);
nav_time   = GPSTime2Date(navData(:,19),navData(:,17), GPSNUMBOOL);
datediff = abs(nav_time-epoch_time);
% rownums = (datediff==min(datediff) & ismember(navData(:,1),PRN));
PRN_idx = ismember(navData(:,1),PRN);
datediff(~PRN_idx) = 1e6;
rownums = (datediff==min(datediff) & ismember(navData(:,1),PRN));

ephemData = navData(rownums,:);

