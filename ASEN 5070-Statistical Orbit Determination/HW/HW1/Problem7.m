% function [X0 Y0 Xp0 Yp0 g] = Problem7()

% Chapter 1, problem 1, Stat OD
% This is basically the equation given in solutions, because I couldn't
% figure it out. 

function Xn=hw1_7
Xn      = [1.5;10;2.2;.5;.3];
j       = zeros(5,1);
deltaj  = zeros(5,5);
rho     = zeros(5,1);
O       = [7;8.00390597;8.94427191;9.801147892;10.630145813];
error   = 1;

while error>10^-9
    for t=0:4
        rho(t+1)=sqrt((Xn(1) + Xn(3)*t -1)^2 + (Xn(2) + Xn(4)*t -Xn(5)*t^2/2 -1)^2);
        deltaj(t+1,1) = (Xn(1) +Xn(3)*t -1)/rho(t+1);
        deltaj(t+1,2) = (Xn(2) +Xn(4)*t - .5*Xn(5)*t^2 -1)/rho(t+1);
        deltaj(t+1,3) =  deltaj(t+1,1)*t;
        deltaj(t+1,4) =  deltaj(t+1,2)*t;
        deltaj(t+1,5) = -deltaj(t+1,2)*t^2/2;
    end

    rho
    deltaj
    j=O-rho;
    XN1 = Xn+inv(deltaj)*j;
    error = norm(XN1 - Xn);
    Xn = XN1;

end








