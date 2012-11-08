function [X_batch,X_kalman,X_Joseph,X_Potter,X_True] = FindStateP4(e)


%% a - Derive P1
% syms err
% H       = [1 2*err;1 3*err];
% R       = eye(2,2);
% P0bar   = (1/err^2)*eye(2,2);
% P1      = inv( transpose(H)*inv(R) * H + inv(P0bar));
% matlabFunction(P1,'file','Find_P1.m');
%% b - Estimate State
Xbar    = [4 2]';
X       = [3 1]'; % True state. Use as observations  
X_True  = X;

%------------------------------------------------------
% Batch Processor
% Xhat    = (H'*inv(R)*H + inv(P1))\(H'*inv(R)*
H       = [1 2*e;1 3*e];
R       = eye(2,2);
P0bar   = (1/e^2)*eye(2,2);
Y       = X;
X_batch = inv(H'*inv(R)*H + inv(P0bar))*(H'*inv(R)*Y + inv(P0bar)*Xbar);


%------------------------------------------------------
% Kalman Filter
H1          = H(1,:);
H2          = H(2,:);
R1          = 1; R2 = R1;
P1bar       = P0bar;  %??????????????????
K1          = P1bar*H1'*inv(H1*P1bar*H1' + R1);
P1_1        = Find_P1(e);
X1_Kalman   = Xbar(1) + K1(1)*(Y(1) - H1*Xbar(1));

% Now again for P2
K2          = P1_1*H2'*inv(H2*P1_1*H2' + R2);
P2_Kalman   = (eye(2) - K2*H2)*P1_1;
X2_Kalman   = Xbar(2) + K2(2)*(Y(2) - H2*Xbar(2));

X_kalman    = [X1_Kalman;X2_Kalman];


%------------------------------------------------------
% Joseph Formulation
K1          = P1bar*H1'*inv(H1*P1bar*H1' + R1);
P1          = (eye(2) - K1*H1)*P1bar*(eye(2)-K1*H1)' + K1*R1*K1';
X1_Joseph   = Xbar(1) + K1*(Y(1) - H1*Xbar(1));

P2bar       = P1;
K2          = P1*H2'*inv(H2*P1*H2' + R2);
P2          = (eye(2) - K2*H2)*P2bar*(eye(2)-K2*H2)' + K2*R2*K2';
X2_Joseph   = Xbar(2) + K2*(Y(2) - H2*Xbar(2));

X_Joseph    = [X1_Joseph;X2_Joseph];


%------------------------------------------------------
% Potter Algorithm
W1bar       = (chol(P0bar))';

Ftilde      = W1bar'*H1';
alpha       = inv(Ftilde'*Ftilde + R1);
gamma       = 1./(1 + sqrt(R1*alpha));
K1           = alpha*W1bar*Ftilde;
% xhat    = xbar + K*([v1;v2]);
W1          = W1bar - gamma*K1*Ftilde';
P1          = W1*W1';
X1_Potter   = Xbar1 + K1*(Y(1) - H1*Xbar(1));

% Now for obs 2
W2bar       = chol(P1)';
Ftilde      = W2bar'*H2';
alpha       = inv(Ftilde'*Ftilde + R2);
gamma       = 1./(1 + sqrt(R2*alpha));
K2          = alpha*W2bar*Ftilde;
% xhat    = xbar + K*([v1;v2]);
W2          = W2bar - gamma*K*Ftilde';
P2_Potter   = W2*W2';
X2_Potter   = Xbar2 + K2*(Y(2) - H2*Xbar(2));


X_Potter    = [X1_Potter;X2_Potter];



