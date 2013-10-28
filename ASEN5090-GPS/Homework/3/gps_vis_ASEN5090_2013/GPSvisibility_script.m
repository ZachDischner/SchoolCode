%==========================================================================
% Simplified driver for GPS Visbility Codes
%
% Based on GPSVisibility_GUI by Ben K. Bradley and calling functions used
% in that GUI.
% P. Axelrad 9/12
%
% Uses:
% generate_GPSyuma_name
% download_GPSyuma
% utc2gpsvec
%==========================================================================
%==========================================================================

clear, close all

% Enter the time of interest
UTC = [2012 9 20 0 0 0];
GPSvec = utc2gpsvec(UTC);


% Construct YUMA almanac file name since this is default setting
navfilename = generate_GPSyuma_name(GPSvec);
[navfilename,statusflag] = download_GPSyuma(GPSvec);


% Grab duration and time step =============================================
durationhrs = 24;

dt_sec = 600;


% Input Antenna Location ==================================================
% North Pole
% latgd = 90;  % latitude, deg
% lon   = 0;  % longitude, deg
% alt   = 0;  % altitude, m

% Equator
% latgd = 0;  % latitude, deg
% lon   = 0;  % longitude, deg
% alt   = 0;  % altitude, m


latgd = 40.0;  % latitude, deg
lon   = -105.0;  % longitude, deg
alt   = 1631.0;  % altitude, m
% % Make antenna pointed straight up
ant_enu = [0 0 1];

% Set minimum and maximum mask angles
mask_min = 10;  % deg
mask_max = 90; % deg

[time_wntow,GPSdata] = ASEN5090_GPSvis(navfilename, 1, GPSvec,...
    durationhrs, dt_sec, latgd, lon, alt,...
    mask_min, mask_max, mask_min, ant_enu, 0, []);
hrofweek = time_wntow(:,2)/3600;

% Plot results ============================================================
% =========================================================================


% Number of Satellites Visible ---------------------------
[ax2] = plot_GPSnumsats(hrofweek,GPSdata.ant_numsats);

title(ax2, 'Number of visible satellites')

% Topocentric: AzEl Plot --------------------------------------------------
[rows,cols] = size(GPSdata.topo_el);

az_vec     = reshape(GPSdata.topo_az,rows*cols,1);
el_vec     = reshape(GPSdata.topo_el,rows*cols,1);
GPSdata.prn = repmat([1:32],rows,1);


% plot the sats visible to antenna only
%fig3 = figure; ax3 = axes;
figure
%plotAzEl(GPSdata.topo_az',GPSdata.topo_el',GPSdata.prn')
plotAzEl(GPSdata.topo_az',GPSdata.topo_el',zeros(rows,cols))

figure
%plotAzEl(GPSdata.topo_az',GPSdata.topo_el',GPSdata.prn')
plotAzEl(GPSdata.topo_az(GPSdata.ant_mask)',GPSdata.topo_el(GPSdata.ant_mask)',zeros(rows,cols))
title('Az and El of just visible satellites')


% Antenna-centric: Elevation Plot -----------------------------------------

time_mat = repmat(hrofweek,1,cols);
time_vec = reshape(time_mat,rows*cols,1);

fig4 = figure; ax4 = axes;
plot(ax4,time_vec,el_vec,'ob','markerfacecolor','b','markersize',4);

ylabel('Elevation (deg)');
xlabel('Time (hr)');

grid(ax4,'on');


title(ax4,{'Elevation Angle';'of Satellites Seen by Antenna'});



