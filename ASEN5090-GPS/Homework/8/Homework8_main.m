%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        Homework8_main.m
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
% type('answers.txt')

%% Setup Work Space
clc;clear all;close all
tic;
screen_size = get(0,'ScreenSize');
sw = screen_size(3);    % Screen Width
sh = screen_size(4);    % Screen Height
addpath HW8_files
soln_format = '|  %2.0f | %15.3f | %7.3f | %12.3f | %15.3f | %9.3f | %3.2f | %3.2f  \t\n';
% h = waitbar(0,'GO PARTY, ILL STAY HERE WORKING!!!');
wb_tot = 1024+1238;
tot_iters = 0;


%% Setup Problem
%------Define Navigation and Observation File
RINEX_FILES = {'onsa2640.onehour','joze2640.onehour'};
nav_msg  = 'brdc2640.12n';

%------Define Orbit determination parameters
params.mu = 3.986005e14;        % Gravitational param [m^3/s^2]
params.we = 7.2921151467e-5;    % Earth's rotation rate [rad/s]

%------Define speed of light
params.c = 299792458; % [m/s]

params.options=optimset('Display','off','TolFun',1e-10,'TolX',1e-10);


%------Define Zenith Correction for each site
Tzenith = [2.3858, 2.4086]; %[m]


%------Read navigation message content
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

%% Calculate  Pre-Fit Residuals
for file_idx=1:length(RINEX_FILES)
    rinex_file = RINEX_FILES{file_idx};
    fprintf('Processing Rinex file %s\n', rinex_file)
    
    %% Read Observation Files, Sort Data
    %------Read a-priori receiver position from header of RINEX obs file
    [ fid, rec_xyz, observables ] = read_rinex_header( rinex_file );
    [rec_lla] = ecef2lla(rec_xyz');
    
    %------Read Observation file
    obs_data    = read_rinex_obs3(rinex_file);
    cols        = obs_data.col; % Structure of column index descriptions
    %------Make nice column addressing variables (P1_col, P2_col ... etc)
    fields = fieldnames(cols);
    for kk=1:length(fields)
        eval(cell2mat([fields(kk),'_col = cols.',fields(kk),';']));
    end
    
    %------Only look at first epoch
    [GPSSecAryUn,secs_idx] = unique(obs_data.data(:,TOW_col));
    t1 = secs_idx(1); t2 = secs_idx(2)-1;
    secs_idx = secs_idx(1);
    GPSSecAryUn = GPSSecAryUn(t1:t2);
    PRNS = obs_data.data(t1:t2 , PRN_col);
    GPSWeekAry = obs_data.data(t1:t2 , WEEK_col);
    GPSSecAry = obs_data.data(t1:t2 , TOW_col);
    %     GPS_Secs = obs_data.data(rows,SOW_col);
    %     GPS_Weeks = obs_data.data(rows,Week_col);
    
    %% Fetch/Compute Pseudorange Values
    
    %------Allocate
    rho_obs = zeros(1, length(PRNS));
    rho_model = rho_obs;
    el = rho_obs;
    res = rho_obs;
    ionFree = rho_obs;
    sat_prn = rho_obs;
    iter = 0;
    for sec_idx = secs_idx(1)'
        
        %% Calculate "Modeled" Pseudorange
        %------Setup Range Finding
        GPS_SOW = GPSSecAry(sec_idx);
        GPS_Week = GPSWeekAry(sec_idx);
        params.Secs = GPS_SOW; %(GPS_Secs(1))  % Seconds used to calculate seconds since epoch
        
        %------Iterate over epoch satellite data
        Nsats = sum(GPSSecAry == GPSSecAry(sec_idx));
        for sat = 1:Nsats
            data_idx = sec_idx+sat-1;
%             waitbar((data_idx+tot_iters)/wb_tot,h)
            iter = iter + 1;
            PRN = obs_data.data(data_idx,PRN_col);
            
            %------Calculate Geometric Range
            [R, rel_dt, satXYZ] = getSatGeomRange(rec_xyz', GPS_Week, GPS_SOW, PRN, nav_data, params);
            
            %------Check Elevation Angle
            [az, el(iter), r] = ecef2azelrange(satXYZ', rec_xyz);
            if el(iter) < 10
                iter=iter-1;
                tot_iters = tot_iters+1;
                continue
            end
            rel_corr = rel_dt*params.c;
            
            %------Get clock correction
            sat_clk_t_corr = getSatClockCorrection(GPS_Week, GPS_SOW, PRN, nav_data);
            
            %------Get Satellite Correction
            sat_corr = sat_clk_t_corr*params.c;
            
            %------Get Tropospheric Correction
            Tropo_corr = getTropoCorrection(Tzenith(file_idx), el(iter) );
            rho_model(iter) = R - sat_corr + Tropo_corr + rel_corr;
            
            %% Fetch Observed Pseudoranges
            %------Get Observed pseudorange
            P2 = obs_data.data(data_idx,P2_col);
            
            if sum(strcmp(observables,'P1'))>0
                %------Retrieve P1 as 2nd freq pseudorange
                PorC1 = obs_data.data(data_idx,P1_col);
%                 rho_obs(iter) = P1;                
            else
                %------Retrieve C1 as 2nd freq pseudorange
                PorC1 = obs_data.data(data_idx,C1_col);
%                 rho_obs(iter) = C1;
            end
            rho_obs(iter) = 2.5457*PorC1-1.5457*P2;
            
            sat_prn(iter) = PRN;
            
%             if iter < secs_idx(2) && file_idx == 1
%                 if sat == 1
%                     fprintf('|_PRN_|___geomRange_____|___rel___|____satClk____|_______P3________|__Prefit___|___el___|__Trop___|\n')
%                 end
%                 fprintf(1,soln_format,PRN,...
%                     R,rel_corr,sat_corr,ionFree(iter),...
%                     ionFree(iter)-rho_model(iter),...
%                     el(iter), Tropo_corr)
%             end
            
        end     % End Satellite Iteration
    end     % End time iteration
    
%     ecef2enuv()
    
    %------Remove data with '0' observation
    zs = rho_obs == 0;
    rho_obs(zs) = []; rho_model(zs) = [];
    %% Plot
    
    
    
%     res = rho_obs-rho_model;
% 
%     
%     obs{file_idx} = rho_obs; %#ok<*SAGROW>
%     model{file_idx} = rho_model;
%     elevation{file_idx} = el;
%     prefit_res{file_idx} = res()';
%     prns = unique(sat_prn);
%     prns(prns==0)=[];
%     
%     %------plot residules
%     figure
%     colors = jet(length(prns));
%     for ii=1:length(prns)
%         rows = sat_prn == prns(ii);
%         sat_res{file_idx, prns(ii)} = res(rows);
%         plot(GPSSecAry(rows),res(rows),'.','color',colors(ii,:),'Markersize',15)
%         leg{ii} = ['PRN: ' , num2str(prns(ii))];
%         hold on
%     end
%     xl = ['Seconds since epoch of week ',num2str(GPS_Week)];
%     xlabel(xl);ylabel('Residual Error [m]')
%     title(rinex_file)
%     legend(leg);
%     
%     %------Plot Histogram
%     figure
%     datan = (res-mean(res))/std(res);
%     hist(datan)
%     histfit(datan)
%     tot_iters = tot_iters + iter;
%     title(rinex_file)
%     
%     %% Ion Free plot for Onsa
%     if strcmp(rinex_file,'onsa2640.onehour')
%         %------Remove data with '0' observation
%         ionFree(zs)=[];
%         %% Plot
%         res = ionFree-rho_model;
%         prns = unique(sat_prn);
%         prns(prns==0)=[];
%         
%         %------plot residules
%         figure
%         colors = jet(length(prns));
%         for ii=1:length(prns)
%             rows = sat_prn == prns(ii);
%             sat_res{file_idx, prns(ii)} = res(rows);
%             plot(linspace(GPSSecAry(1),GPSSecAry(end),length(res(rows))),res(rows),'.','color',colors(ii,:),'Markersize',15)
%             leg{ii} = ['PRN: ' , num2str(prns(ii))];
%             hold on
%         end
%         xl = ['Seconds since epoch of week ',num2str(GPS_Week)];
%         xlabel(xl);ylabel('Residual Error [m]')
%         title(['Ion-Free residules for ', rinex_file])
%         legend(leg);
%         
%         figure
%         for ii=1:length(prns)
%             rows = sat_prn == prns(ii);
%             sat_res{file_idx, prns(ii)} = res(rows);
%             plot(el(rows),res(rows),'.','color',colors(ii,:),'Markersize',15)
%             leg{ii} = ['PRN: ' , num2str(prns(ii))];
%             hold on
%         end
%         xlabel('Elevation Angle [Degrees]');ylabel('Residule [m]')
%         title('Ion Free Residule vs Elevation Angle for Onsa')
%         legend(leg)
%         
%         %------Plot Histogram
%         figure
%         datan = (res-mean(res))/std(res);
%         hist(datan)
%         histfit(datan)
%         title(['Ion Free normalized hist: ', rinex_file])
%         
%     end
%     
%     clear leg rows outliers sat_prn
    
end     % End Rinex File Iteration


%% Clean, Reformat
figure_awesome


%% SUPPORTING FUNCTION - Homework8_main.m
type('Homework8_main.m')

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
