%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        getStationPosition.m
% Author      : Zach Dischner
% Date        : 11/18/2013
% Description : Calculate a reciever station position given input RINEX and
%                   nav files
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
% Inputs        :
%
% Outputs       :
%
% TODOS         :
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

function [recXYZ_ECI, del_x] = getStationPosition(nav_data, obs_data, ...
                            rec_xyz, observables, el_mask, Tzenith, params)

%% Read Observation Files, Sort Data

% [rec_lla] = ecef2lla(rec_xyz');

%------Form least squares like variables
xr = rec_xyz(1); yr = rec_xyz(2); zr = rec_xyz(3);
X = [xr, yr, zr, 1];

%------Sort Obs variables
cols        = obs_data.col; % Structure of column index descriptions

%------Make nice column addressing variables (P1_col, P2_col ... etc)
fields = fieldnames(cols);
for kk=1:length(fields)
    eval(cell2mat([fields(kk),'_col = cols.',fields(kk),';']));
end

%------Only look at first epoch
[GPSSecAryUn,secs_idx] = unique(obs_data.data(:,TOW_col));
PRNS = obs_data.data(: , PRN_col);
GPSWeekAry = obs_data.data(: , WEEK_col);
GPSSecAry = obs_data.data(: , TOW_col);

%% Fetch/Compute Pseudorange Values

%------Allocate
rho_obs = zeros(1, length(PRNS));
rho_model = rho_obs;
el = rho_obs;
res = rho_obs;
ionFree = rho_obs;
sat_prn = rho_obs;
epoch_iter = 0;
for sec_idx = secs_idx'
    epoch_iter = epoch_iter + 1;
    
    %% Calculate "Modeled" Pseudorange
    %------Setup Range Finding
    GPS_SOW = GPSSecAry(sec_idx);
    GPS_Week = GPSWeekAry(sec_idx);
    params.Secs = GPS_SOW; %(GPS_Secs(1))  % Seconds used to calculate seconds since epoch
    
    %------Iterate over epoch satellite data
    Nsats = sum(GPSSecAry == GPSSecAry(sec_idx));
    iter = 0; clear H y
    for sat = 1:Nsats
        data_idx = sec_idx+sat-1;
        iter = iter + 1;
        PRN = obs_data.data(data_idx,PRN_col);
        
        %------Calculate Geometric Range
        [R, rel_dt, satXYZ] = getSatGeomRange(rec_xyz', GPS_Week, GPS_SOW, PRN, nav_data, params);
        xs = satXYZ(1); ys = satXYZ(2); zs = satXYZ(3);
        
        %------Check Elevation Angle
        [az, el(iter), r] = ecef2azelrange(satXYZ', rec_xyz);
        if el(iter) < el_mask
            iter=iter-1;
            continue
        end
        rel_corr = rel_dt*params.c;
        
        %------Get clock correction
        sat_clk_t_corr = getSatClockCorrection(GPS_Week, GPS_SOW, PRN, nav_data);
        
        %------Get Satellite Correction
        sat_corr = sat_clk_t_corr*params.c;
        
        %------Get Tropospheric Correction
        Tropo_corr = getTropoCorrection(Tzenith, el(iter) );
        rho_model(iter) = R - sat_corr + Tropo_corr + rel_corr;
        
        %% Fetch Observed Pseudoranges
        P2 = obs_data.data(data_idx,P2_col);
        
        if sum(strcmp(observables,'P1'))>0
            %------Retrieve P1 as 2nd freq pseudorange
            PorC1 = obs_data.data(data_idx,P1_col);
        else
            %------Retrieve C1 as 2nd freq pseudorange
            PorC1 = obs_data.data(data_idx,C1_col);
        end
        
        %------Use ONLY ionosphere free. Uncomment linesto use what data we have
        if P2 == 0 || PorC1 == 0
            iter=iter-1;
            continue
        else
            rho_obs(iter) = 2.5457*PorC1-1.5457*P2;
        end
        
        %------Populate Least Squares Variables
        H(iter,:) = [(xr-xs)/R, (yr-ys)/R, (zr-zs)/R, 1]; %??? What about clock???
        % same as (rec_xyz' - satXYZ)/R
        y(iter) = rho_obs(iter) - rho_model(iter);
        
        sat_prn(iter) = PRN;
        
%         if epoch_iter == 1 && file_idx == 1 && debug == 1
%             if sat == 1
%                 fprintf('|_PRN_|________P3_______|___Geom Range___|____satClk____|__rel__|__Tropo___|___prefit___|\n')
%             end
%             fprintf(1,soln_format,PRN,...
%                 rho_obs(iter), R,sat_corr,...
%                 rel_corr,Tropo_corr,y(iter))
%         end
        
    end     % End Satellite Iteration
    %------Obtain least squares correction
    del_x(epoch_iter,:) = transpose((H'*H)\H'*y');
    P = inv(H'*H).*(0.5)^2;
    
    %------Rotate
    rec_xyz_adj = repmat(rec_xyz', epoch_iter, 1) + del_x(:,1:3);
    rec_lla_adj = ecef2lla(rec_xyz_adj);
    lat_station =  rec_lla_adj(epoch_iter,1);
    lon_station = rec_lla_adj(epoch_iter,2);
    P(1:3,1:3) = rotateCovariance(P(1:3,1:3), lat_station,lon_station);
    sigmas(epoch_iter,:) = sqrt(diag(P))';
    
    %         key_epoch = 2;
    %         if epoch_iter == key_epoch
    %             fprintf('\n\n----------------------------%s-------------------------',rinex_file)
    %             rec_xyz_adj = repmat(rec_xyz', epoch_iter, 1) - del_x(:,1:3);
    %             rec_lla_adj = ecef2lla(rec_xyz_adj);
    %             [de, dn, du] = ecef2enuv(del_x(key_epoch,1), ...
    %                                      del_x(key_epoch,2), ...
    %                                      del_x(key_epoch,3), ...
    %                                      rec_lla_adj(key_epoch,1),rec_lla_adj(key_epoch,2));
    %            fprintf('\nEast/Sig (m) %3.3f %3.3f\n',de,sigmas(epoch_iter,1))
    %            fprintf('North/Sig (m) %3.3f %3.3f\n',dn,sigmas(epoch_iter,2))
    %            fprintf('Up/Sig (m) %3.3f %3.3f\n',du,sigmas(epoch_iter,3))
    %            fprintf('Clk/Sig (m) %5.3f %3.3f\n',del_x(key_epoch,4),sigmas(epoch_iter,4))
    %         end
end     % End time iteration


%------Finish up and return
recXYZ_ECI = rec_xyz_adj;






