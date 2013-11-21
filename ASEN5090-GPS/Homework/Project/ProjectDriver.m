%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        ProjectDriver.m
% Authors     : Zach Dischner
%               Pierce Martin
%               Greg Nelson
%               Andrew Haynes
% Date        : 11/17/2013
% Description : Matlab driver script for ASEN 5090 Group Project
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
tic;
addpath SupportingFiles ObservationFiles
% soln_format = '| %3.0f  |  %12.8f  |  %12.8f  |  %12.4f  |  \t\n';
debug = 0;


%% Setup Problem
%------Define Navigation and Observation File
RINEX_FILES = {'8-11-11_nya12200.11o'};
nav_msg  = '8-11-11_nyal2200.11n';

%------Define Orbit determination parameters
params.mu = 3.986005e14;        % Gravitational param [m^3/s^2]
params.we = 7.2921151467e-5;    % Earth's rotation rate [rad/s]

%------Define speed of light
params.c = 299792458; % [m/s]

params.options=optimset('Display','off','TolFun',1e-20,'TolX',1e-20);

%------Elevation Mask
el_mask = 10; %degrees


%------Define Zenith Correction for each site
Tzenith = [2.3, 2.3]; %[m]

%------Define convergence criterion
conv_tol = 1e-6;
noconv_count = 100;

%------Read in nav message here (for speed purposes)
nav_data = read_GPSbroadcast(nav_msg); % Returns [n x 25] matrix of sat orbit information

%% Find Station Positions
for file_idx=1:length(RINEX_FILES)
%     rinex_file = RINEX_FILES{file_idx};
    count = 0;
    delta = 1000;
    
    %------Read a-priori receiver position from header of RINEX obs file
    [ fid, rec_xyz_priori, observables ] = read_rinex_header( RINEX_FILES{file_idx} );
    
    %------Read Obs data (for speed purposes, do it here)
    obs_data    = read_rinex_obs3(RINEX_FILES{file_idx});
    
    fprintf('\n\n----------------------------%s-------------------------\n',RINEX_FILES{file_idx})
    
    while delta > conv_tol
        %% Iterate on Receiver Position Until Convergance
        [recXYZ_ECI, del_x, model] = getStationPosition(nav_data, obs_data, ...
                                        rec_xyz_priori, observables, ...
                                        el_mask, Tzenith(file_idx), params);
        recLLA = ecef2lla(recXYZ_ECI);
        
        if count == 0
            fprintf('\n\n|_iter_|______lat_______|_______lon______|_____height_____|\n')
        end
        fprintf(1,soln_format, count, recLLA(1),recLLA(2),recLLA(3))
        
        if count > noconv_count
            fprintf(2,'RECEIVER POSITION NOT CONVERGING!!')
            break
        end
        
        %------Compute current delta
        delta(count+1) = norm(recXYZ_ECI' - rec_xyz_priori);
        rec_xyz_priori = recXYZ_ECI';
        
        count = count + 1;
    end
    lat(file_idx) = recLLA(1);
    lon(file_idx) = recLLA(2);
end  % End Rinex File iteration




% %% Clean, Reformat

% %% SUPPORTING FUNCTION - Homework8_main.m
% type('Homework9_main.m')
% 
% %% SUPPORTING FUNCTION - getStationPosition.m
% type('getStationPosition.m')
% 
% %% SUPPORTING FUNCTION - getSatGeomRange.m
% type('getSatGeomRange.m')
% 
% %% SUPPORTING FUNCTION - date2GPSTime.m
% type('date2GPSTime.m')
% 
% %% SUPPORTING FUNCTION - findNearestEphem.m
% type('findNearestEphem.m')
% 
% %% SUPPORTING FUNCTION - calculateSatellitePosition.m
% type('calculateSatellitePosition.m')
% 
% %% SUPPORTING FUNCTION - findFirstEpoch.m
% type('findFirstEpoch.m')
% 
% %% SUPPORTING FUNCTION - getSatClockCorrection.m
% type('getSatClockCorrection.m')
% 
% %% SUPPORTING FUNCTION - date2GPSTime.m
% type('GPSTime2Date.m')
% 
% %% SUPPORTING FUNCTION - getTropoCorrection.m
% type('getTropoCorrection.m')

fprintf('\nSim took %3.1f seconds to run\n',toc)
