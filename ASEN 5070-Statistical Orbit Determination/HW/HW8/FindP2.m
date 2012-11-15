function [P2_True,P2_Kalman,P2_Joseph,P2_Potter,P2_Batch]=FindP2(e,tr)
% Compute P2 (or the trace therof) for various different algorithms. 

% Input is error value, and 'tr', indicator to return the trace of each 
%       calculation


if exist('tr','var')
    tr=1;
else 
    tr = 0;
end
    
% Return the trace!
% eq 4.7.20
% z1 = |1 3||x1| + |v1|
% z2 = |1 1||x2| + |v2|

% Looks like    y=Hx+e
P1bar   = (1/e^2)*eye(2,2);
H1      = [1 e];
H2      = [1 1];
R       = 1;
%------------------------------------------------------------
%% a - True Solution from 4.7.24

B   = 1 - 2*e + 2*e^2*(2+e^2);
P2_True  = 1/B*[ 1+2*e^2     -(1+e); ...
                  -(1+e)       2+e^2 ];
if tr              
    P2_True     = trace(P2_True);
end
%------------------------------------------------------------

%% b - P2 Using Conventional Kalman Filter
K1 = P1bar*H1'*inv(H1*P1bar*H1' + R);
P1 = (eye(2) - K1*H1)*P1bar;

% Now again for P2
K2          = P1*H2'*inv(H2*P1*H2' + R);
P2_Kalman   = (eye(2) - K2*H2)*P1;

% P2_Kalman =  1/(1-2*e)*[   -1  1; ...
%                          1   -1];
if tr
    P2_Kalman    = trace(P2_Kalman);
end
% matlabFunction(P2_Kalman,'file','P2_Kalman.m');   

%------------------------------------------------------------

%% c - P2 Using Joseph Formulation
% Branching off of the Kalman
K1          = P1bar*H1'*inv(H1*P1bar*H1' + R);
P1          = (eye(2) - K1*H1)*P1bar*(eye(2)-K1*H1)' + K1*R*K1';
K2          = P1*H2'*inv(H2*P1*H2' + R);

P2_Joseph   = (eye(2) - K2*H2)*P1*(eye(2)-K2*H2)' + K2*R*K2';
if tr
    P2_Joseph     = trace(P2_Joseph);         
end

%------------------------------------------------------------

%% d - P2 Using Plotter Algorithm (5.7.17) 

W1bar   = (chol(P1bar))';

Ftilde  = W1bar'*H1';
alpha   = inv(Ftilde'*Ftilde + R);
gamma   = 1./(1 + sqrt(R*alpha));
K       = alpha*W1bar*Ftilde;
% xhat    = xbar + K*([v1;v2]);
W1      = W1bar - gamma*K*Ftilde';
P1      = W1*W1';

% Now for obs 2
W2bar   = chol(P1)';
Ftilde  = W2bar'*H2';
alpha   = inv(Ftilde'*Ftilde + R);
gamma   = 1./(1 + sqrt(R*alpha));
K       = alpha*W2bar*Ftilde;
% xhat    = xbar + K*([v1;v2]);
W2      = W2bar - gamma*K*Ftilde';

P2_Potter = W2*W2';
if tr
    P2_Potter     = trace(P2_Potter);
end

%------------------------------------------------------------

%% e - P2 Using Batch
H       = [1 e;1 1];        %H       = [1 e; 1 1];   
R       = eye(2,2);            %   eye(2,2);      % E(ee')
P2bar   = (1/e)^2*eye(2,2);
Delta   = P2bar\eye(2,2);
Delta   = Delta + transpose(H)*inv(R)*H;

P2_Batch= Delta\eye(2,2);
if tr
    P2_Batch= trace(P2_Batch);
end

          