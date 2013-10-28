%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        Homework5_main.m
% Author      : Zach Dischner
% Date        : 10/10/2013
% Description : Matlat script for all calculations required for
%               ASEN 5090 Homework 5
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
% screen_size = get(0,'ScreenSize');
% sw = screen_size(3);    % Screen Width
% sh = screen_size(4);    % Screen Height
% figColor = [0.99 0.99 0.98];
addpath HW5_files

%% Setup Problem
%------Define Navigation and Observation File
nav_msg  = 'brdc2640.12n';
obs_file = 'test.12o';

%------Read navigation message content
nav_ephem = read_GPSbroadcast(nav_msg); % Returns [n x 25] matrix of sat orbit information
%                  col1: prn, PRN number of satellite
%                  col2: M0, mean anomaly at reference time, rad
%                  col3: delta_n, mean motion difference from computed value, rad/s
%                  col4: ecc, eccentricity of orbit
%                  col5: sqrt_a, square root of semi-major axis, m^0.5
%                  col6: Loa, longitude of ascending node of orbit plane at weekly epoch, rad
%                  col7: incl, inclination angle at reference time, rad
%                  col8: perigee, argument of perigee, rad
%                  col9: ra_rate, rate of change of right ascension, rad/s
%                 col10: i_rate, rate of change of inclination angle, rad/s
%                 col11: Cuc, amplitude of the cosine harmonic correction term to the argument of latitude
%                 col12: Cus, amplitude of the sine harmonic correction term to the argument of latitude
%                 col13: Crc, amplitude of the cosine harmonic correction term to the orbit radius
%                 col14: Crs, amplitude of the sine harmonic correction term to the orbit radius
%                 col15: Cic, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col16: Cis, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col17: Toe, reference time ephemeris (seconds into GPS week)
%                 col18: IODE, issue of data (ephemeris) 
%                 col19: GPS_week, GPS Week Number (to go with Toe)
%                 col20: Toc, time of clock
%                 col21: Af0, satellite clock bias (sec)
%                 col22: Af1, satellite clock drift (sec/sec)
%                 col23: Af2, satellite clock drift rate (sec/sec/sec)
%                 col24: blank (zero)
%                 col25: health, satellite health (0=good and usable)

%% Get GPS Time from UTC time
%------Convert time of interest into GPS week and seconds of week
% Using standard for GPS Epoch: http://tycho.usno.navy.mil/gps_week.html
% Verified with online calendar: http://adn.agi.com/GNSSWeb/
emph_date      = 'September 20 2012 02:11:00';
[GPS_W, GPS_SOW] = date2GPSTime(emph_date);
fprintf('\n*^*^*^*^*^*^*^*^*^*^*^*^*^*^* HW5 Step 3 *^*^*^*^*^*^*^*^*^*^*^*^*^*^*\n')
fprintf('GPS Week and seconds of week: [ %d , %d ]\n',GPS_W, GPS_SOW)



%% Calculate Satellite Position
%------Define Orbit determination parameters
params.mu = 3.986005e14;    % Gravitational param [m^3/s^2]
params.we = 7.2921151467e-5;    % Earth's rotation rate [rad/s]
params.Secs = GPS_SOW;      % Seconds used to calculate seconds since epoch

%------Define Satellite PRN Ranges
PRNs = [15, 17, 12];
fprintf('\n\n*^*^*^*^*^*^*^*^*^*^*^*^*^*^* HW5 Step 6 *^*^*^*^*^*^*^*^*^*^*^*^*^*^*\n\n')
for PRN = PRNs
    %------Fetch corresponding ephemeris data
    ephem = findNearestEphem(PRN,GPS_W,GPS_SOW, nav_ephem);
    [x,y,z] = calculateSatellitePosition(ephem,params);
    fprintf('PRN: %d ===> [x,y,z]: [%10.3f , %10.3f , %10.3f]m \n',PRN,x,y,z) 
    [h,pos,v,r,s,a]=broadcast2xva(nav_ephem,[GPS_W,GPS_SOW],PRN);
    fprintf('   Ben B Ref [x,y,z]: [%10.3f , %10.3f , %10.3f]m \n\n',pos)
    
end
fprintf('\n\n')
    

%% SUPPORTING FUNCTION - date2GPSTime.m
type('date2GPSTime.m')

%% SUPPORTING FUNCTION - findNearestEphem.m
type('findNearestEphem.m')

%% SUPPORTING FUNCTION - calculateSatellitePosition.m
type('calculateSatellitePosition.m')
