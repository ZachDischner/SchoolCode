% plot_X_star.m
% Create a plot of the satellite trajectory
% Author: Pierce Martin
% Inputs:
%   X_star: Full 18 element state at each observation time (Each column
%   should correspond to the state at each obs time, ie size 18x385)
%   RE: Radius of earth
%   station_id: Array of station ids at each obs time
function [] = plot_X_star(X_star,RE,station_id)

% Create figure to plot into
figure, hold on

% Load the basic MATLAB earth topographic data
load('topo.mat','topo','topomap1');

% Clear the axis
cla reset

% Create the earth surface.
[x,y,z] = sphere(50);
x = RE*x*0.99;
y = RE*y*0.99;
z = RE*z*0.99;
props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;
load topo
zlimits = [min(topo(:)) max(topo(:))];
demcmap(zlimits);

shading('interp')
surface(x,y,z,props);


% Add lights.
light('position',[-1 0 1]);
light('position',[-1.5 0.5 -0.5], 'color', [.6 .2 .2]);

% Plot the satellite trajectory
hold on
plot3(X_star(1,:),X_star(2,:),X_star(3,:),'.r');

% Plot the stations
plot3(X_star(10,end),X_star(11,end),X_star(12,end),'.g','MarkerSize',20,'LineWidth',4)
plot3(X_star(13,end),X_star(14,end),X_star(15,end),'.g','MarkerSize',20,'LineWidth',4)
plot3(X_star(16,end),X_star(17,end),X_star(18,end),'.g','MarkerSize',20,'LineWidth',4)

% Plot lines from station to satellite (use every 5th obs)
for i = 1:5:size(X_star,2)
    
    % Determine which station the observation is from
    switch station_id(i)
        
        case 101
            xs = X_star(10,end);
            ys = X_star(11,end);
            zs = X_star(12,end);
            
        case 337
            xs = X_star(13,end);
            ys = X_star(14,end);
            zs = X_star(15,end);
            
        case 394
            xs = X_star(16,end);
            ys = X_star(17,end);
            zs = X_star(18,end);
           
            
    end
    
    % Create the line
    line([X_star(1,i) xs],[X_star(2,i) ys],[X_star(3,i) zs]);
    
%     if mod(i,10) == 0
%         error_ellipse(P(1:3,1:3,i),[X_star(1,i),X_star(2,i),X_star(3,i)]);
%     end

end

% Set the view.
axis square off

end