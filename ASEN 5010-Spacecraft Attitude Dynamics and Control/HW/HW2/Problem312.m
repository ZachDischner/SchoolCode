function Q = Problem312()
clc;clear all;close all
%% Initial Conditions
Q_init  = deg2rad([40 30 80]);
% Time space 
t=linspace(0,60,1000);

%% Perform the integration
Q = zeros(length(t),length(Q_init));
Q(1,:)= Q_init;
for ii=2:length(t)
   % Find the Derivative
   Qprime = IntYPR(t(ii),Q(ii-1,:));
   
   % General linear integration:
   %    x_(n+1) = x_(n) + x'*delta_t
   Q(ii,:) = Q(ii-1,:) + Qprime'*(t(ii)-t(ii-1));
   
end

 [t,Q]=ode45(@IntYPR,t,Q_init');
 
domod=0;

figure;hold on
xlabel('Time [s]');ylabel('Angle ({$^\circ$})')
title('Time Evolution of $\psi , \theta , \phi$ Angles')
plot(t,mod(rad2deg(Q(:,1)),360*domod),'b');
plot(t,mod(rad2deg(Q(:,2)),360*domod),'r');
plot(t,mod(rad2deg(Q(:,3)),360*domod),'g');
legend('{$\Psi $}','{$\theta $}','{$\phi $}','location','best')
end



function Qprime = IntYPR(t,Q)
    % Calculate [B(th)] Matrix
    w=Calcw(t);
    B = BmatEuler321(Q);
    Qprime = B*w;
end

function w = Calcw(t)
w       =  [sin(0.1*t) ; 
            0.01 ; 
            cos(0.1*t)].*deg2rad(20);
end