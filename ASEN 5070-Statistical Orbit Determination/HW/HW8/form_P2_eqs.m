function form_P2_eqs()
% eq 4.7.20
% z1 = |1 3||x1| + |v1|
% z2 = |1 1||x2| + |v2|

% Looks like    y=Hx+e

%------------------------------------------------------------
%% a - True Solution from 4.7.24
syms B e v1 v2
B   = 1 - 2*e + 2*e^2*(2+e^2);
P2_true  = 1\B*[ 1+2*e^2     -(1+e); ...
                  -(1+e)       2+e^2 ]; 
matlabFunction(P2_true,'file','P2_true.m');

%------------------------------------------------------------

%% b - P2 Using Conventional Kalman Filter
P2_Kalman =  1/(2-e)*[   -1  1; ...
                         1   -1]; 
matlabFunction(P2_Kalman,'file','P2_Kalman.m');   

%------------------------------------------------------------

%% c - P2 Using Joseph Formulation
P2_Joseph = [ 1+2*e    -(1+3*e);...
             -(1+3*e) 2+e];
matlabFunction(P2_Joseph,'file','P2_Joseph.m');

%------------------------------------------------------------

%% d - P2 Using Plotter Algorithm (5.7.17) 
H       = [1 e; 1 1];
Wbar    = inv((1/e)^2*eye(size(H)));
R       = [v1 0;0 v2];

Ftilde  = transpose(Wbar)*transpose(H);
alpha   = inv(transpose(Ftilde)*Ftilde + R);
gamma   = 1/(1 + sqrt(R*alpha));
K       = alpha*Wbar*Ftilde;
xhat    = xbar + K*([v1;v2]);
W       = Wbar - gamma*K*transpose(Ftile);

P2_Potter = inv(W);



%------------------------------------------------------------

%% e - P2 Using Batch
H       = [1 e; 1 1];   
R       = [v1 0;0 v2];
P1bar   = (1/e)^2*eye(size(H));
Delta   = inv(P1bar);
Delta   = Delta + transpose(H)*inv(R)*H;

eval(['P2_Batch = @(e,v1,v2) ',inv(Delta),';']); 
P2_Batch = [ (v1 + e^2*v2 + e^2*v1*v2)/(2*e^2*v1 - 2*e + e^2*v2 + e^4*v2 + e^2 + e^4*v1*v2 + 1),          -(v1 + e*v2)/(2*e^2*v1 - 2*e + e^2*v2 + e^4*v2 + e^2 + e^4*v1*v2 + 1);
              -(v1 + e*v2)/(2*e^2*v1 - 2*e + e^2*v2 + e^4*v2 + e^2 + e^4*v1*v2 + 1), (v1*v2*e^2 + v1 + v2)/(2*e^2*v1 - 2*e + e^2*v2 + e^4*v2 + e^2 + e^4*v1*v2 + 1)];
matlabFunction(P2_Batch,'file','P2_Batch.m');          
          
          