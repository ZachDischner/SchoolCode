function [A Htilde] = FindA_Htilde()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Zach Dischner-10/21/2012
% 
% FindA_Htilde
% 
% Purpose: Find A and Htilde matrices
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% Use Jacobian

%% Define the State
syms x y z xdot ydot zdot uE J2 Cd Xsite1  Ysite1  Zsite1 Xsite2 Ysite2 Zsite2  Xsite3 Ysite3 Zsite3 theta theta_dot
syms R_E r Area m rho_a Va va 

X       = [x ; y ; z ; xdot ; ydot ; zdot ; uE ; J2 ; Cd ; Xsite1 ; Ysite1 ; Zsite1; Xsite2; Ysite2; Zsite2 ; Xsite3; Ysite3; Zsite3];

%% Define F vector

% [xdot ydot zdot] due to gravity and J2
r       = sqrt(x^2+y^2+z^2);
F_U     = [...
              -uE/r^3*x*(1-3/2*J2*(R_E/r)^2*(5*(z/r)^2-1));
              -uE/r^3*y*(1-3/2*J2*(R_E/r)^2*(5*(z/r)^2-1));
              -uE/r^3*z*(1-3/2*J2*(R_E/r)^2*(5*(z/r)^2-3));
          ];
      
% [xdot ydot zdot] due to atmospheric drag
Va      = [
    xdot + theta_dot*y;
    ydot - theta_dot*x;
    zdot
    ];
va      = sqrt((xdot + theta_dot*y)^2 + (ydot - theta_dot*x)^2 + zdot^2);

syms rho0 r0 H
rho_a   = rho0.*exp((r - r0)./H);

F_Drag = -0.5 .* Cd .* (Area./m) .* rho_a .* va .* Va;

% Assemble
F_a = F_U - F_Drag;

% X' = F*X
F =[xdot ; ydot ; zdot ; F_a ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];
      

%% Find A Matrix
clear r; syms r
A = jacobian(F,X);
% Do this outside of the function?
A=simplify(A);
A=subs(A,sqrt(x^2+y^2+z^2),r);
A=subs(A,(x^2+y^2+z^2),r^2);
A = simplify(A);


%% Find Htilde Matrix (for a generic station position)
syms Xsite Ysite Zsite
syms theta  % = theta_dot * time
% dcm = [cos -sin 0;sin cos 0; 0 0 1];
% rho = sqrt((x-Xsite*cos(theta) - Ysite*sin(theta))^2 + (y-Ysite*cos(theta) + Xsite*sin(theta))^2 + (z-Zsite)^2);
% The Answer
rho = sqrt(x^2+y^2+z^2+Xsite^2+Ysite^2+Zsite^2+2*(x*Xsite+y*Ysite)*cos(theta)+2*(x*Ysite-y*Xsite)*sin(theta)-2*z*Zsite);
% rho = @(x,y,z,Xsite,Ysite,Zsite,theta) sqrt(x^2+y^2+z^2+Xsite^2+Ysite^2+Zsite^2+2*(x*Xsite+y*Ysite)*cos(theta)+2*(x*Ysite-y*Xsite)*sin(theta)-2*z*Zsite)

rhodot = (x*xdot + y*ydot + z*zdot - (xdot*Xsite + ydot*Ysite)*cos(theta) + theta_dot*(x*Xsite + y*Ysite)*sin(theta)...
            +(xdot*Ysite - ydot*Xsite)*sin(theta) + theta_dot*(x*Ysite - y*Xsite)*cos(theta) - zdot*Zsite)...
                                                            /rho;
obs = [rho;rhodot];

% For just single station
Single_Station_X       = [x ; y ; z ; xdot ; ydot ; zdot ; uE ; J2 ; Cd ; Xsite ; Ysite ; Zsite];

% Add padding depending on the station. 
% If station 1, add 6 columns of padding
% If station 2, add 3 columns, move 3 endmost columns after that, then add
% 3 columns
% IF station 3, put 6 columns of zeros between the 9th and (going to 15th) columms. 

Htilde = jacobian(obs,Single_Station_X);

                             



















end
