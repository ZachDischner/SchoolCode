%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        GPSTime2Date.m
% Author      : Zach Dischner
% Date        : 10/24/2013
% Description : Convert a date type object into [GPS_Weeks, GPS_SOW] time
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
% Inputs        : [GPS_Weeks, GPS_SOW]-weeks and seconds of week
% 
% Outputs       : utcDate - Satellite PRN number
% 
% TODOS         : Vectorize!
%                 Build in mod options
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function utcDate  = GPSTime2Date(GPS_Weeks, GPS_SOW)
% gps_week_start = 'January 6 1980 00:00:00';
gps_weeks_start = 723186;   % datenum(gps_week_start) Save time
%------GPS date in numerical form, since Matlab's 'epoch'
GPS_Num = (GPS_Weeks+GPS_SOW/7/24/3600)*7 + gps_weeks_start;
utcDate = datestr(datevec(GPS_Num));