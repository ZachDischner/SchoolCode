function y = intEOM2(t,X)
m = 1;
c = 0.1;
k1=1;
k3 = 0.1 ;
P = 5;

x = X(1);
xdot = X(2);

xr = sin(0.5*t);
xrdot = 0.5*cos(0.5*t);
xrddot = -0.25*sin(0.5*t);

u = m*xrddot - k3*(x-xr)^3 + k3*x^3 + c*xdot + k1*xr - P*(xdot-xrdot);

xddot = 1/m * (u - c*xdot - k1*x - k3*x^3);

y = [xdot;xddot];
