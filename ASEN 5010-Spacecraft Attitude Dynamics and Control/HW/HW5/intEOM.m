function y = intEOM(t,X)
m = 1;
c = 0.1;
k1=1;
k3 = 0.1 ;
P = 5;

x = X(1);
xdot = X(2);

xr = 0;
xrdot = 0;
xrddot = 0;

u = m*xrddot - k3*(x-xr)^3 + k3*x^3 + c*xdot + k1*xr - P*(xdot-xrdot);

xddot = 1/m * (u - c*xdot - k1*x - k3*x^3);

y = [xdot;xddot];

