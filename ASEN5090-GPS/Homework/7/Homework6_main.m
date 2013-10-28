%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        Homework6_main.m
% Author      : Zach Dischner
% Date        : 10/24/2013
% Description : Matlab script for all calculations required for
%               ASEN 5090 Homework 6
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
% figColor = [0.99 0.99 0.98];
addpath HW6_files
soln_format = '|  %2.0f | %15.3f | %7.3f | %12.3f | %15.3f | %7.3f \t\n';

%% Setup Problem
%------Define Navigation and Observation File
nav_msg  = 'brdc2640.12n';
obs_file = 'darw264x.12o';
fprintf('1) Navigation File: %s\n2)Observation File: %s\n\n',nav_msg,obs_file);

%------Define Orbit determination parameters
params.mu = 3.986005e14;        % Gravitational param [m^3/s^2]
params.we = 7.2921151467e-5;    % Earth's rotation rate [rad/s]

%------Define speed of light
params.c = 299792458; % [m/s] 

%% Read Files
%------Read navigation message content
fprintf('3) Read Navigation File\n\n')
nav_data = read_GPSbroadcast(nav_msg); % Returns [n x 25] matrix of sat orbit information
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

%------Read a-priori receiver position from header of RINEX obs file
fprintf('4) Get a-priori from RINEX file\n\n')
[ fid, rec_xyz, observables ] = read_rinex_header( obs_file );

%------Read Observation file
obs_data = read_rinex_obs3(obs_file);
Week_col = 1;
SOW_col = 2;    % Simple indicator for clarification
PRN_col = 3;    % Simple indicator for clarification
C1_col = 6;
rows = find(obs_data.data(:,SOW_col)==min(obs_data.data(:,SOW_col)));
PRNS = obs_data.data(rows,PRN_col);
GPS_Secs = obs_data.data(rows,SOW_col);
GPS_Weeks = obs_data.data(rows,Week_col);



%% Calculate Geometric Range for First Epoch Satellites
fprintf('5) Get ephemeris data for first epoch in rinex file\n\n')
[epochData,rows] = findNearestEphem(PRNS,GPS_Weeks(1),GPS_Secs(1),nav_data);

fprintf(['6)For all the PRNs in the first epoch, make (and call)', ...
         'a function \n\tthat calculates the geomet- ric range (use instructions',...
                    '\n\tat the end of this assignment). Since your broadcast ',...
                    '\n\tephemeris has the information needed, calculate the ',...
                    '\n\trelativity correction.\n\n'])
type('getSatGeomRange')
fprintf('7) Write a function that calculates satellite clock correction\n\n')
type('getSatClockCorrection.m')
fprintf('8) Access values for C1\n\t[>>C1(ii) = obs_data.data(ii,C1_col);]\n\n')
fprintf('9) Output values in readable format\n')

%------Allocate
Tt = zeros(length(rows),1);
R=Tt; sat_clk_t_corr=Tt; satcorr=Tt; rel_corr=Tt; C1=Tt;
fprintf('|_PRN_|___geomRange_____|___rel___|____satClk____|________C1_______|__C1-R+satcorr\n')
for ii = 1:length(rows)
    %------Setup Range Finding
    GPS_SOW = epochData(ii,17);
    GPS_Week = GPS_Weeks(1);
    params.Secs = GPS_Secs(1);  % Seconds used to calculate seconds since epoch
    
    %------Calculate Geometric Range
    [R(ii), rel_dt] = getSatGeomRange(rec_xyz', GPS_Week, GPS_Secs(1), PRNS(ii), nav_data, params);
    rel_corr(ii) = rel_dt*params.c;
    %------Get clock correction
    sat_clk_t_corr(ii) = getSatClockCorrection(GPS_Week, GPS_Secs(1), PRNS(ii), nav_data); 
    
    %------Get Satellite Correction
    satcorr(ii) = sat_clk_t_corr(ii)*params.c;
 
    %------Retrieve C1
    C1(ii) = obs_data.data(ii,C1_col);

        
    %------Output Answers yo!
    fprintf(1,soln_format,PRNS(ii),...
        R(ii),rel_corr(ii),satcorr(ii),C1(ii),C1(ii)-R(ii)+satcorr(ii))
end



%% SUPPORTING FUNCTION - date2GPSTime.m
type('date2GPSTime.m')

%% SUPPORTING FUNCTION - findNearestEphem.m
type('findNearestEphem.m')

%% SUPPORTING FUNCTION - calculateSatellitePosition.m
type('calculateSatellitePosition.m')

%% SUPPORTING FUNCTION - findFirstEpoch.m
type('findFirstEpoch.m')

%% SUPPORTING FUNCTION - date2GPSTime.m
type('date2GPSTime.m')
