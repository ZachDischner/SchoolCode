 function [GPSdata] = gps_masks(GPSdata, sz, topomaskmin, topomaskmax, antmask, ant_enu)

%% Compute Antenna Az, El vectors
azAnt = wrapTo360( atan2d( ant_enu(2), ant_enu(1) ) );
elAnt = asind(ant_enu(3)/norm(ant_enu));

%% Trig Calculations for antenna
SazAnt = sind(azAnt);
CazAnt = cosd(azAnt);

SelAnt = sind(90 - elAnt);
CelAnt = cosd(90 - elAnt);

%% Trig Calculations for satellite
Saz2 = sind(GPSdata.topo_az);
Caz2 = cosd(GPSdata.topo_az);

Sel2 = sind(90 - GPSdata.topo_el);
Cel2 = cosd(90 - GPSdata.topo_el);

%% Angle between sat and ant, used by dot product
alpha = acosd((CazAnt.*Caz2 + SazAnt.*Saz2).*SelAnt.*Sel2 + CelAnt.*Cel2);

%% Satellite visibility matrix
topo_numsats = zeros(sz,32);
ant_numsats = zeros(sz,32);


%% Create Topographical Mask
topo=GPSdata.topo_el;
topo(topo ~= topo) = 0;  % Take care of nans
topo(topo < topomaskmin) = 0;
topo(topo > topomaskmatopo) = 0;
topo(topo ~= 0) = 1;
GPSdata.topo_numsats = sum(topo,2);

%% Create Antenna Mask
c1 = GPSdata.topo_el > topomaskmin;
c2 = GPSdata.topo_el < topomaskmax;
c3 = alpha < 90 - antmask;              % Adjust for spherical Trig

% Apply conditions -- satellite is visible if all are true
topo_numsats(c1 & c2) = 1;
ant_numsats( c1 & c2 & c3) = 1;

% Sum to get number of sats visible
GPSdata.ant_numsats  = sum(ant_numsats,2);

end %function