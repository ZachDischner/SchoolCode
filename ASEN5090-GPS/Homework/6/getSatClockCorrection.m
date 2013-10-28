function [tcorr] = getSatClockCorrection(GPS_Weeks, GPS_SOW, PRN, nav_data) 
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        getSatClockCorrection.m
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
% Inputs        : PRN - Satellite PRN number
%                 GPSWeeks - GPS week number (modded or no?)
%                 GPSSOW - GPS Seconds of week
%                 navData - A full array of all emphimeris data, fetched
%                            from navigation file
% Outputs       : t_corr - Satellite clock correction
% 
% History       Oct 24 2013 - First Rev
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%------Get ephemeris dataset
[eph_data,tmp] = findNearestEphem(PRN, GPS_Weeks, GPS_SOW, nav_data);

%------Define readibility indices
Af0_col = 21;   %Af0, satellite clock bias (sec)
Af1_col = 22;   %Af1, satellite clock drift (sec/sec)
Af2_col = 23;   %Af2, satellite clock drift rate (sec/sec/sec)
SOW_col = 17;   %Toe, reference time ephemeris (seconds into GPS week)

%------Fetch Correction Constants
Af0 = eph_data(Af0_col); 
Af1 = eph_data(Af1_col); 
Af2 = eph_data(Af2_col);

t_eph = eph_data(SOW_col);
dt = GPS_SOW - t_eph;

%------Calculate clock correction
tcorr = Af0 + Af1*(dt) + Af2*(dt)^2;

end %function



