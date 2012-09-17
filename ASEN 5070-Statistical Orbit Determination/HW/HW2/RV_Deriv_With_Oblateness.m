% Derivative function, designed to return the derivative of the Cartesian
%   R and V observations, incliuding the effect of the Earth's Oblateness
%
% Written by Zach Dischner on 8/30/2012
%
% inputs:
%         time: time...
%         RV  : A column vector containing R and V measurements as
%                       X
%                       Y
%                       Z
%                       X'
%                       Y'
%                       Z'
%
%
% outputs: (All in ECI I think)
%         RVdot: The time derivative of R and V vectors. (I.E. [V;A])
%                       X'
%                       Y'
%                       Z'
%                       X''
%                       Y''
%                       Z''
%
%
%


function [RVdot] = RV_Deriv_With_Oblateness(time,RV)

u   = 398600.4;     % km^3/s^2
J2  = 0.00108248;   % []
R_E = 6378.145;     % km    Radius of earth

if size(RV) ~= [6 1]
    % For easy computation of VA later on down the line
    RVdot = [RV(:,4:end),[(-u.*RV(:,1)./sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2).^3)' ; (-u.*RV(:,2)./sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2).^3)' ; (-u.*RV(:,3)./sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2).^3)']'];

else
    x = RV(1);
    y = RV(2);
    z = RV(3);
    r = sqrt(x^2 + y^2 + z^2);
    
    RVdot = [...
                                             RV(4:end);        
                               -u/r^3*x*(1-3/2*J2*(R_E/r)^2*(5*(z/r)^2-1));
                               -u/r^3*y*(1-3/2*J2*(R_E/r)^2*(5*(z/r)^2-1));
                               -u/r^3*z*(1-3/2*J2*(R_E/r)^2*(5*(z/r)^2-3));

        ];
    
end


end


