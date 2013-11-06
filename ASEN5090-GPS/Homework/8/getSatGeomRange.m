%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        getSatGeomRange.m
% Author      : Zach Dischner
% Date        : 10/22/2013
% Description : calculate satellite position from GPS ephemeris dataset
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
% Inputs        : rStation - GPS Rx [x,y,z] coords in ECEF meters
%                 GPS_Weeks - GPS Week time
%                 GPS_SOW - Seconds into week
%                 PRN - Satellite PRN
%                 nav_data - nx25 array of sat data from broadcase
%                   ephemeris 
%                 params - structure containing keplarian specs and extra
%                           calculations for sat position
% 
% Outputs       : R - Geometric Range value, in meters
%                 rel_dt - clock offset due to relativity (in seconds)
%                 rSat - 3d ECI coordinates of satellite [x,y,z]
% 
% History       October 11 2013 - First Rev
%               October 24 2013 - Reformatted output to [R,tk]
%                               - Added check for time field in params,
%                                 other that that in the ephemeris data
%               October 30 2013 - Reformatted output to include XYZ sat
%                                 position
%               
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [R, rel_dt, rSat] = getSatGeomRange(rStation, GPS_Weeks, GPS_SOW, PRN, nav_data, params)

%------Find Nearest Ephemeris
[epochData,rows] = findNearestEphem(PRN, GPS_Weeks, GPS_SOW, nav_data);
SOW_col = 20;
% Single Row in this case

%------Get Sat Position from Ephemeris data
[rSat,rel_dt] = calculateSatellitePosition(epochData, params);


%------Set up convergence limits
R = 0;
conv_limit = 1e-15;
max_iters = 200;
iter = 1;

%------Iterate and converge on Geometric Range
while(1)
    %------Calculate Geometric Range
    Rtmp = norm( rSat - rStation );
    
    %------Check for Convergence
    if(abs(Rtmp - R) < conv_limit)
        break
    end
    
    %------Assign new Range Value now that criterion are passed
    R = Rtmp;
    
    %------Check for iteration limit
    if(iter > max_iters)
        fprintf(2,'YO BRO!! Range Calculation not converging!\n')
        break
    end
    
    %------Increase iteration count
    iter = iter + 1;
    
    %------Calculate 'Tt', time of transmission
    dt = R/params.c;
%     Tr = epochData(SOW_col);
    Tr = GPS_SOW;
    Tt = Tr - dt;
    
    %------Recalculate Satellite position
    params.Secs = Tt;     
    % to use new time value
    [rSat,rel_dt,params] = calculateSatellitePosition(epochData,params);
    
    %------Rotate Sat position at time Tr (account for earth's rotation)
    phi = params.we*dt;
    rSat = transpose(rot3(phi)*rSat');
    
end
rmfield(params,'E_guess');

