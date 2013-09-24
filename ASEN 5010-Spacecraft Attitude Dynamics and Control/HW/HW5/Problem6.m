clc;clear all;close all;
m = 1;
c = 0.1;
k1=1;
k3 = 0.1 ;

x0 = 2;
xdot0 = 0;


t = linspace(0,25);
tol=1e-15;
[T,X] = ode45(@intEOM,t,[x0;xdot0],['reltol',tol,'abstol',tol]);

figure(1)
xr1 = zeros(length(t),1);
plot(t,xr1,'k.');grid on;
hold on
plot(t,X(:,1),'k-')
xlabel('Time (s)');
ylabel('x (m)');
title('Performance on Reference 1, P = 5');
legend('Reference','Control Performance');


figure(2)
[T,X2] = ode45(@intEOM2,t,[x0;xdot0],[]);
xr2 = sin(0.5*t);
plot(t,xr2,'k.'); grid on;
hold on
plot(t,X2(:,1),'k-')
xlabel('Time (s)');
ylabel('x (m)');
title('Performance on Reference 2, P = 5');
legend('Reference','Control Performance');

