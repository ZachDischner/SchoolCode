function [X_batch,X_Kalman,X_Joseph,X_Potter,X_True] = FindStateP4(e)
% Return state calculations with various algorithms for homework 8,
% problem 4. 


%% b - Estimate State
Xbar    = [4 2]';
X       = [3 1]'; % True state. Use as observations  
X_True  = X';

%------------------------------------------------------
% Batch Processor
% Xhat    = (H'*inv(R)*H + inv(P1))\(H'*inv(R)*
H       = [1 2*e;1 3*e];
R       = eye(2,2);
P0bar   = (1/e^2)*eye(2,2);
Y       = H*X;
X_batch = inv(H'*inv(R)*H + inv(P0bar))*(H'*inv(R)*Y + inv(P0bar)*Xbar);
X_batch = X_batch';


%------------------------------------------------------
% Kalman Filter - Think of the two observations as occuring at different
% times. With the first becoming a-priori for the next one
H1          = H(1,:);
H2          = H(2,:);
R1          = 1; R2 = R1;
P1bar       = P0bar; 

K1          = P1bar*H1'*inv(H1*P1bar*H1' + R1);
P1_1        = (eye(2) - K1*H1)*P1bar;
X1_Kalman   = Xbar + K1*(Y(1) - H1*Xbar);

% Now again for number 2. The X1_Kalman state is now the a-priori for the
% new state
K2          = P1_1*H2'*inv(H2*P1_1*H2' + R2);
P2_Kalman   = (eye(2) - K2*H2)*P1_1;
X2_Kalman   = X1_Kalman + K2*(Y(2) - H2*X1_Kalman);

X_Kalman    = X2_Kalman';


%------------------------------------------------------
% Joseph Formulation
K1          = P1bar*H1'*inv(H1*P1bar*H1' + R1);
P1          = (eye(2) - K1*H1)*P1bar*(eye(2)-K1*H1)' + K1*R1*K1';
X1_Joseph   = Xbar + K1*(Y(1) - H1*Xbar);

P2bar       = P1;
K2          = P1*H2'*inv(H2*P1*H2' + R2);
P2          = (eye(2) - K2*H2)*P2bar*(eye(2)-K2*H2)' + K2*R2*K2';
X2_Joseph   = X1_Joseph + K2*(Y(2) - H2*X1_Joseph);

X_Joseph    = X2_Joseph';


%------------------------------------------------------
% Potter Algorithm
W1bar       = (chol(P0bar))';

Ftilde      = W1bar'*H1';
alpha       = inv(Ftilde'*Ftilde + R1);
gamma       = 1./(1 + sqrt(R1*alpha));
K1          = alpha*W1bar*Ftilde;
% xhat    = xbar + K*([v1;v2]);
W1          = W1bar - gamma*K1*Ftilde';
P1          = W1*W1';
X1_Potter   = Xbar + K1*(Y(1) - H1*Xbar);

% Now for obs 2
W2bar       = chol(P1)';
Ftilde      = W2bar'*H2';
alpha       = inv(Ftilde'*Ftilde + R2);
gamma       = 1./(1 + sqrt(R2*alpha));
K2          = alpha*W2bar*Ftilde;
% xhat    = xbar + K*([v1;v2]);
W2          = W2bar - gamma*K2*Ftilde';
P2_Potter   = W2*W2';
X2_Potter   = X1_Potter + K2*(Y(2) - H2*X1_Potter);


X_Potter    = X2_Potter';



