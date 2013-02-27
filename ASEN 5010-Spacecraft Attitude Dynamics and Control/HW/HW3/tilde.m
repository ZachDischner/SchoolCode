%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        tilde.m
% Author      : Zach Dischner
% Date        : 2/25/2012
% Description : Perform tilde operation on 3 element vector
%
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

function xtilde = tilde(x)
    xtilde = [  0   -x(3)   x(2);
               x(3)   0    -x(1);
              -x(2)  x(1)    0  ];

end