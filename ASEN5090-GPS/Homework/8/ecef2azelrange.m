function [az,el,range] = ecef2azelrange(r_sat,r_site)

%==========================================================================
%==========================================================================
% [az,el,range] = ecef2azelrange(r_sat,r_site,latgd,lon)
%
% Calculates the azimuth, elevation, and range of a satellite with respect
% to an observation site. At a specific instance in time.
%
%
% Author: Pierce Martin
% Date: 10/26/2013
%
% INPUT:         Description                                         Units
%
%  r_sat      - position of satellite in ECEF frame               [x y z] m
%  r_site     - position of observing site in ECEF frame          [x y z] m
%
%
% OUTPUT:       
%    
%  az         - azimuth (degrees clockwise from North)          [0,360] deg
%  el         - elevation (degrees up from horizon)            [-90,90] deg
%  range      - distance from observation site to satellite              m                             
%
%
% Coupling:
%
%  none
%
%
%==========================================================================
%==========================================================================

% Calculate current latitude and longitude of the reciever
lla = ecef2lla(r_site');
latgd = lla(1);
lon = lla(2);

% Calculate UP unit vector in ECEF coordinates
e_up_ECEF = [cosd(latgd)*cosd(lon); cosd(latgd)*sind(lon); sind(latgd)];

% Calculate line of site vector from the site to the satellite in ECEF
e_LOS_ECEF = (r_sat - r_site)/norm(r_sat - r_site);

% Calculate rotation matrix from ECEF to ENU
ecef2enu = [-sind(lon)              cosd(lon)               0;
            -sind(latgd)*cosd(lon)  -sind(latgd)*sind(lon)  cosd(latgd);
            cosd(latgd)*cosd(lon)   cosd(latgd)*sind(lon)   sind(latgd)];
        
% Calculate line of site vector in ENU coords
e_LOS_ENU = ecef2enu*e_LOS_ECEF;
           
% Compute elevation, azimuth and range
az = atan2d(e_LOS_ENU(1),e_LOS_ENU(2));
el = asind(dot(e_up_ECEF,e_LOS_ECEF));
range = norm(r_sat - r_site);