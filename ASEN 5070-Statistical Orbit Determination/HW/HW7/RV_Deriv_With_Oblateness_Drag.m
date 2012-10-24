% Derivative function, designed to return the derivative of the Cartesian
%   R and V observations, incliuding the effect of the Earth's Oblateness
%
% Written by Zach Dischner on 8/30/2012
%
% inputs:
%         time: time...
%         State  : A column vector containing R and V measurements as
%                       X
%                       Y
%                       Z
%                       X'
%                       Y'
%                       Z'
%                       u_e
%                       J2
%                       Cd
%                       Xbar_sitei   (1 2 or 3) based on station coordinates
%                       Phi
%
%
% outputs: (All in ECI I think)
%         StateDot: The time derivative of R and V vectors. (I.E. [V;A])
%                       X'
%                       Y'
%                       Z'
%                       X''
%                       Y''
%                       Z''
%
%
%


function [RVdot] = RV_Deriv_With_Oblateness_Drag(time,State)

u   = 398600.4;     % km^3/s^2
J2  = 0.00108248;   % []
R_E = 6378.145;    % km    Radius of earth
% 
% if size(State) ~= [6 1]
%     % For easy computation of VA later on down the line
%     RVdot = [State(:,4:end),[(-u.*State(:,1)./sqrt(State(:,1).^2 + State(:,2).^2 + State(:,3).^2).^3)' ; (-u.*State(:,2)./sqrt(State(:,1).^2 + State(:,2).^2 + State(:,3).^2).^3)' ; (-u.*State(:,3)./sqrt(State(:,1).^2 + State(:,2).^2 + State(:,3).^2).^3)']'];
% 
% else
    
    % Accl due to gravity and Oblateness
    x  = State(1);
    y  = State(2);
    z  = State(3);
    vx = State(4);
    vy = State(5);
    vz = State(6);
    r  = sqrt(x^2 + y^2 + z^2);
    
    
    A_U =[...
                -u/r^3*x*(1-3/2*J2*(R_E/r)^2*(5*(z/r)^2-1));
                -u/r^3*y*(1-3/2*J2*(R_E/r)^2*(5*(z/r)^2-1));
                -u/r^3*z*(1-3/2*J2*(R_E/r)^2*(5*(z/r)^2-3));
         ];
   
   
    % Accl due to Drag 
    CD      = 2.0;
    Area    = 3.6;                  % m^2
    m       = 1350;                 % kg
    rho0    = 4e-13;                % kg/m^3
    r0      = 7298.145;             % km
    H       = 200;                  % km
    TH_dot  = 7.29211585530066e-5;  % rad/s
    
    rho_a   = rho0.*exp((r - r0)./H);
    Va      = [
                vx + TH_dot*y;
                vy - TH_dot*x;
                     vz
              ];
          
    va      = norm(Va);
    
    A_Drag = -0.5 .* CD .* (Area./m) .* rho_a .* va .* Va; 
   
    A = A_U - A_Drag;
    
    RVdot = [...
                                             State(4:end);
                                                 A


        ];
    
end




