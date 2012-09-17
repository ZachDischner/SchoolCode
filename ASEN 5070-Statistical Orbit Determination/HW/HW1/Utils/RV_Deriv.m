% Derivative function, designed to return the derivative of the Cartesian
%   R and V observations
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


function [RVdot] = RV_Deriv(time,RV)

u   = 398600.5;     %km^3/s^2

if size(RV) ~= [6 1]
    
    RVdot = [RV(:,4:end),[(-u.*RV(:,1)./sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2).^3)' ; (-u.*RV(:,2)./sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2).^3)' ; (-u.*RV(:,3)./sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2).^3)']'];
%                                 
%                         RV(:,4:end),
%            [(-u.*RV(:,1)./sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2).^3)' ;
%            (-u.*RV(:,2)./sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2).^3)' ; 
%            (-u.*RV(:,3)./sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2).^3)']'...
%            ];
% %         %                             -RV(1:3).*u./(norm(RV(1:3)).^3);
%         [(-u.*RV(:,1)./sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2).^3)',
%          (-u.*RV(:,2)./sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2).^3)',
%          (-u.*RV(:,3)./sqrt(RV(:,1).^2 + RV(:,2).^2 + RV(:,3).^2).^3')]'
%         ];
    
else
    RVdot = [...
                                RV(4:end);
        %                             -RV(1:3).*u./(norm(RV(1:3)).^3);
        -u.*RV(1)./sqrt(RV(1).^2 + RV(2).^2 + RV(3).^2).^3;
        -u.*RV(2)./sqrt(RV(1).^2 + RV(2).^2 + RV(3).^2).^3;
        -u.*RV(3)./sqrt(RV(1).^2 + RV(2).^2 + RV(3).^2).^3;
        ];
    
end


end


