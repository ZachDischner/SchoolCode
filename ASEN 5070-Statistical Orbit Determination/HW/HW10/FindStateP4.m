% function [X_Givens,X_True] = FindStateExtras(epsilon)
% Return state calculations with various algorithms for homework 10
% e is \delta in homework description
clc;clear all;close all;format long G
epsilon = 1e-8;
tic
%% b - Estimate State
Xbar    = [4 7]';
xbar    = Xbar;   %????
X       = [3 1]'; % True state. Use as observations 
X_obs   = [3 2]';
y       = X_obs;   % Book Notation

X_True  = X';
H       = [1 epsilon;1 1];
h       = H;       % Book Notation
R       = eye(2,2);
P0bar   = (1/epsilon^2)*eye(2,2);

% Givens Algorithm
%------------------------------------------------------
% Find a-priori stuff for  Givens
% Choleski decompisition
Sbar    = (chol(P0bar))';

% Find Rbar
Rbar    = inv(Sbar);       % Not R

% Find d
dbar       = diag(Rbar).^2;

% Find Ubar
for ii=1:length(Rbar(:,1))
    for jj=1:length(Rbar(1,:))
        Ubar(ii,jj) = Rbar(ii,jj)/Rbar(ii,ii);
    end
end

Dbar = (Rbar*inv(Ubar)).^2;

% Find btilde (= btilde_bar?)
b_bar = Rbar*xbar;
tmp = Dbar.^(-1/2);
tmp(tmp == inf) = 0;
btilde_bar  = tmp*b_bar;
btilde = btilde_bar;


m       = length(X_obs);
n       = length(P0bar);   %  P  --> Sbar --> R --> n

Sum = 0;
% U       = zeros(size(Ubar));
U       = Ubar;
U(1,1)  = 1;
d       = dbar;
for kk = 1:m
   delta_k = 1;
   for ii = 1:n
      d_prime(ii)   = d(ii) + delta_k*(h(kk,ii)^2);
      Cbar          = d(ii)/d_prime(ii);
      Sbar          = delta_k*(h(kk,ii))/d_prime(ii);
      y_prime(kk)   = y(kk) - btilde(ii)*h(kk,ii);
      btilde(ii)    = btilde(ii)*Cbar + y(kk)*Sbar;
      y(kk)         = y_prime(kk);
      delta_k       = delta_k*Cbar;
      d(ii)         = d_prime(ii);
       for jj = ii+1:n
           h_tilde(kk,jj) = h(kk,jj) - U(ii,jj)*h(kk,ii);
           U(ii,jj) = U(ii,jj)*Cbar + h(kk,jj)*Sbar;
           h(kk,jj) = h_tilde(kk,jj);
       end
       
       D(ii,ii) = d(ii);
   end

   e(kk)    =sqrt(delta_k*y(kk));
   Sum = Sum + e(kk)^2;
end


xhat_Givens=U\btilde

P2_Givens = inv(U)*inv(D)*inv(U)'


%% True
Beta = 1-2*epsilon + 2*epsilon^2*(2+epsilon^2);
P2_True = 1/Beta*[ 1+2*epsilon^2 -(1+epsilon); -(1+epsilon) 2+epsilon^2]

toc

% %------------------------------------------------------
% % Batch Processor
% % Xhat    = (H'*inv(R)*H + inv(P1))\(H'*inv(R)*
% Y       = H*X;
% X_batch = inv(H'*inv(R)*H + inv(P0bar))*(H'*inv(R)*Y + inv(P0bar)*Xbar);
% X_batch = X_batch';
% 
% 
% %------------------------------------------------------
% % Kalman Filter - Think of the two observations as occuring at different
% % times. With the first becoming a-priori for the next one
% H1          = H(1,:);
% H2          = H(2,:);
% R1          = 1; R2 = R1;
% P1bar       = P0bar; 
% 
% K1          = P1bar*H1'*inv(H1*P1bar*H1' + R1);
% P1_1        = (eye(2) - K1*H1)*P1bar;
% X1_Kalman   = Xbar + K1*(Y(1) - H1*Xbar);
% 
% % Now again for number 2. The X1_Kalman state is now the a-priori for the
% % new state
% K2          = P1_1*H2'*inv(H2*P1_1*H2' + R2);
% P2_Kalman   = (eye(2) - K2*H2)*P1_1;
% X2_Kalman   = X1_Kalman + K2*(Y(2) - H2*X1_Kalman);
% 
% X_Kalman    = X2_Kalman';
% 
% 
% %------------------------------------------------------
% % Joseph Formulation
% K1          = P1bar*H1'*inv(H1*P1bar*H1' + R1);
% P1          = (eye(2) - K1*H1)*P1bar*(eye(2)-K1*H1)' + K1*R1*K1';
% X1_Joseph   = Xbar + K1*(Y(1) - H1*Xbar);
% 
% P2bar       = P1;
% K2          = P1*H2'*inv(H2*P1*H2' + R2);
% P2          = (eye(2) - K2*H2)*P2bar*(eye(2)-K2*H2)' + K2*R2*K2';
% X2_Joseph   = X1_Joseph + K2*(Y(2) - H2*X1_Joseph);
% 
% X_Joseph    = X2_Joseph';
% 
% 
% %------------------------------------------------------
% % Potter Algorithm
% W1bar       = (chol(P0bar))';
% 
% Ftilde      = W1bar'*H1';
% alpha       = inv(Ftilde'*Ftilde + R1);
% gamma       = 1./(1 + sqrt(R1*alpha));
% K1          = alpha*W1bar*Ftilde;
% % xhat    = xbar + K*([v1;v2]);
% W1          = W1bar - gamma*K1*Ftilde';
% P1          = W1*W1';
% X1_Potter   = Xbar + K1*(Y(1) - H1*Xbar);
% 
% % Now for obs 2
% W2bar       = chol(P1)';
% Ftilde      = W2bar'*H2';
% alpha       = inv(Ftilde'*Ftilde + R2);
% gamma       = 1./(1 + sqrt(R2*alpha));
% K2          = alpha*W2bar*Ftilde;
% % xhat    = xbar + K*([v1;v2]);
% W2          = W2bar - gamma*K2*Ftilde';
% P2_Potter   = W2*W2';
% X2_Potter   = X1_Potter + K2*(Y(2) - H2*X1_Potter);
% 
% 
% X_Potter    = X2_Potter';
















