%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        Homework9_main.m
% Author      : Zach Dischner
% Date        : 11/3/2013
% Description : Matlab script for all calculations required for
%               ASEN 5090 Homework 8
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
%% Questions to Answer (First for ease of grading)
type('answers.txt')

%% Setup Work Space
clc;clear all;close all
tic;
screen_size = get(0,'ScreenSize');
sw = screen_size(3);    % Screen Width
sh = screen_size(4);    % Screen Height
addpath HW9_files
soln_format = '| %3.0f  |  %12.8f  |  %12.8f  |  %12.4f  |  \t\n';
debug = 0;


%% Setup Problem
%------Define Navigation and Observation File
RINEX_FILES = {'xxxx2640.12o.txt','yyyy2640.12o.txt'};
nav_msg  = 'brdc2640.12n';

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
parfor file_idx=1:length(RINEX_FILES)
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
        [recXYZ_ECI, del_x] = getStationPosition(nav_data, obs_data, ...
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
    %% Show on a Map Because I'm a Visual Creature
%     name = 'Zach Dischner Consulting';
%     description = sprintf('GPS Location for file [%s]',rinex_file);
%     webmap();
%     wmmarker(lat(file_idx),lon(file_idx),'featureName',rinex_file,'Description',description,'OverlayName',name)
%     if file_idx == 1
%         wmzoom(13);
%     else
%         wmzoom(15);
%     end
end  % End Rinex File iteration

% webmap();
% description = sprintf('GPS Location for file [%s]',RINEX_FILES{1});
% wmmarker(lat(1),lon(1),'featureName',RINEX_FILES{1},'Description',description,'OverlayName',name)
% description = sprintf('GPS Location for file [%s]',RINEX_FILES{2});
% wmmarker(lat(2),lon(2),'featureName',RINEX_FILES{2},'Description',description,'OverlayName',name)



% %% Clean, Reformat

%% SUPPORTING FUNCTION - Homework8_main.m
type('Homework9_main.m')

%% SUPPORTING FUNCTION - getStationPosition.m
type('getStationPosition.m')

%% SUPPORTING FUNCTION - getSatGeomRange.m
type('getSatGeomRange.m')

%% SUPPORTING FUNCTION - date2GPSTime.m
type('date2GPSTime.m')

%% SUPPORTING FUNCTION - findNearestEphem.m
type('findNearestEphem.m')

%% SUPPORTING FUNCTION - calculateSatellitePosition.m
type('calculateSatellitePosition.m')

%% SUPPORTING FUNCTION - findFirstEpoch.m
type('findFirstEpoch.m')

%% SUPPORTING FUNCTION - getSatClockCorrection.m
type('getSatClockCorrection.m')

%% SUPPORTING FUNCTION - date2GPSTime.m
type('GPSTime2Date.m')

%% SUPPORTING FUNCTION - getTropoCorrection.m
type('getTropoCorrection.m')

fprintf('\nSim took %3.1f seconds to run\n',toc)
