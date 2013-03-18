clc;clear all;close all;

EA  = deg2rad([-30 40 20])';
C   = Euler3132C(EA);

%% a- Find ehat princ rot axis
fprintf('Finding PRV \n\n')
% ehat= C2PRV(C);
%  Problem with this function????
Phi = acos(.5*(C(1,1)+C(2,2)+C(3,3)-1));

ehat=1/(2*sin(Phi)).*[C(2,3)-C(3,2);
                      C(3,1)-C(1,3);
                      C(1,2)-C(2,1)]
if ehat ~= C2PRV(C)/Phi; fprintf('ALERT! Schaubs C2PRV Doesnt Match!!!\n\n'); end;                  
                  
%% b - Phi, Phi'
fprintf('Finding Phi, Phi'' \n\n')
Phi
Phiprime=Phi-2*pi


%% c - Quaternions baby!!!
fprintf('Finding Quaternions \n\n')
b=[ cos(Phi/2);
    ehat(1)*sin(Phi/2);
    ehat(2)*sin(Phi/2);
    ehat(3)*sin(Phi/2)
  ]

if b ~= C2EP(C); fprintf('ALERT! Schaubs C2EP Doesnt Match!!!\n\n'); end;

%% d - CRP
fprintf('Finding CRP \n\n')

CRP = [ b(2:4)]/b(1)

if CRP ~= C2Gibbs(C); fprintf('ALERT! C2gibbs Toolbox Doesnt Match!!!\n\n'); end;  

%% e - MRP
fprintf('Finding CRP \n\n')

MRP = [ b(2:4)]/(1+b(1))

if MRP ~= C2MRP(C); fprintf('ALERT! C2gibbs Toolbox Doesnt Match!!!\n\n'); end;  




