function Prot = rotateCovariance(P,lat,lon)
% input P(1:3,1:3). Just xyz coord vectors

% Calculate rotation matrix from ECEF to ENU
ecef2enu = [-sind(lon)              cosd(lon)               0;
            -sind(lat)*cosd(lon)  -sind(lat)*sind(lon)  cosd(lat);
            cosd(lat)*cosd(lon)   cosd(lat)*sind(lon)   sind(lat)];
        
Prot = ecef2enu*P*ecef2enu';
        
        
end