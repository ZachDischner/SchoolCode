clc;clear all; close all
% Problem 3.14, HW2
% Zach Dischner
ehat=1/sqrt(3)*[1 1 1]';
e1=ehat(1);
e2=ehat(2);
e3=ehat(3);

phi=pi/4;

% fprintf('Shaubs direct PRV2Euler321:\n')
% YPR=PRV2Euler321([phi;ehat])
% fprintf('Doesnt Match!!!\n\n')

fprintf('HUsing Schaubs methods:\n')

sig=1-cos(phi);

C=[ e1^2*sig+cos(phi)       ,   e1*e2*sig+e3*sin(phi)   ,   e1*e3*sig-e2*sin(phi)   ;
    e2*e1*sig-e3*sin(phi)   ,     e2^2*sig+cos(phi)     ,   e2*e3*sig+e1*sin(phi)   ;
    e3*e1*sig+e2*sin(phi)   ,   e3*e2*sig-e1*sin(phi)   ,     e3^2*sig+cos(phi)     ];

YPR2=C2Euler321(C);

fprintf('Y P R angles are:\n\n')
rad2deg(YPR2)

% All By Hand
Y=atan2(C(1,2),C(1,1));
P=-asin(C(1,3));
R=atan2(C(2,3),C(3,3));

fprintf(' Or by doing a hand calculation:\n\n')
rad2deg([Y,P,R])   

fprintf('\nBoth match up\n\n')
