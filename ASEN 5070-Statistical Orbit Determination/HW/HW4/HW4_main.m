clc;clear all;close all; format compact
% Go back and comment the shit out of this!!!

%% 1
% 1a
syms x y k

k = 1/int(int((x^2+y^2),0,2),1,3);
k = eval(k);
disp(['k is:  ',num2str(k)]);

% 1b
clear x; clear y; 
fun     = @(x,y) k*( x.^2 + y.^2);
ymin    = 2;
ymax    = 3;
xmin    = 1; % maybe 0
xmax    = 2;
P_1B    = quad2d(fun,xmin,xmax,ymin,ymax);
% P =int(int(k*(x^2 + y^2),2,3),1,2);
disp(['Probability 1B:  ',num2str(P_1B.*100),'%']);

% 1C

% P =int(int(k*(x^2 + y^2),1,3),1,2);
P_1C    = quad2d(fun,1,3,1,2);
disp(['Probability 1B:  ',num2str(P_1C.*100),'%']);

% 1D
fun     = @(x,y) k*( x.^2 + y.^2);
ymin    = @(x) 4-x;
ymax    = 3;
xmin    = 1; 
xmax    = 2;

P_1D    = quad2d(fun,xmin,xmax,ymin,ymax);
disp(['Probability 1D:  ',num2str(P_1D.*100),'%']);

% 1E
xmin    = 0;
xmax    = 1;
ymin    = @(x) 4-x;
ymax    = @(x) 4-x;
P_1E    = quad2d(fun,xmin,xmax,ymin,ymax);
disp(['Probability 1E:  ',num2str(P_1E.*100),'%     (0 expected, no volume)']);

% 1F
fun_x   = @(x) k*(x.^2 + (3)^2);
Hy      = @(y) k*( 2*y.^2 + 8/3);
xmin    = 0;    xmax    = 1;
y_value = 3;
P_1F    = quad(fun_x,xmin,xmax)/Hy(y_value);
disp(['Probability 1F:  ',num2str(P_1F.*100),'%']);

%1G -Standard Deviation
% ???????
% marginal distribution, how to do this in matlab cleanly?
Gx      = @(x) (2*k*(3*x.^2 + 13))/3;
E_x2    =  quad(@(x) (x.^2).*Gx(x),0,2);
E_x_2   = (quad(@(x) x.*Gx(x),0,2))^2;

disp(['E(x^2):  ',num2str(E_x2)]);
disp(['E(x)^2:  ',num2str(E_x_2)]);
Sigma_x = sqrt(E_x2 - E_x_2);
disp(['sigma_x:  ',num2str(Sigma_x)]);
    % PROBABLY WRONG!!!
% SigmaX  




% 1H
% marginal density function. Wrong???
% Hy = int(k*(x.^2+y.^2),x,0,2);


P_1H    = quad2d(fun,1,2,1,2)./quad(Hy,1,2);
disp(['Probability 1H:  ',num2str(P_1H.*100),'%']);






 
