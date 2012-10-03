function [yprime] = RV_Deriv1b(time,y)

% Simple function to pass to integrator for HW5, problem 1b. 
% yprime = d(y)/dt = 

% where y = [x y x' y' PHI]
% 
% so yprime = [x' y' x'' y'' PHI']
% 
% where x'' = -x/r^3;  y'' = -y/r^3; phi' = A*phi

u = 1;



r = sqrt(y(1).^2 + y(2).^2);
 
A = [ 1     0       0       0
      0     1       0       0
      0     0    -1./r.^3   0;
      0     0       0   -1./r.^3;
    ];

A = [               0                       0                       1        0;
                    0                       0                       0        1;
    -u./(r.^3)+3*u.*y(1).^2./(r.^5)   3*u*y(1)*y(2)/(r^5)           0        0;
            3*u*y(1)*y(2)/(r^5)    -u./(r.^3)+3*u.*y(2).^2./(r.^5)  0        0;
    ];
      

% lumping together RV' = A*y     and PHI' = A*PHI

% PHI = reshape(y(5:20),4,4)
      
yprime = [[y(3) y(4) -y(1)/r^3 -y(2)/r^3 ]'; reshape( A * reshape(y(5:20),4,4), 16,1)];

end