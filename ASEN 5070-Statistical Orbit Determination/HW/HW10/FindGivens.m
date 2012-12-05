function [P_Givens,X_Givens,P_True,Sum_e,Sum_e_NEW] = FindStateExtras(epsilon)
% Return state calculations with various algorithms for homework 10
% e is \delta in homework description
% clc;clear all;close all;format long G
% epsilon = 1e-8;
% tic
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

Sum_e = 0;
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

   e(kk)    =sqrt(delta_k)*y(kk);
   Sum_e = Sum_e + e(kk)^2;
end

xhat_Givens=U\btilde;


P2_Givens = inv(U)*inv(D)*inv(U)';



%% True
Beta = 1-2*epsilon + 2*epsilon^2*(2+epsilon^2);
P2_True = 1/Beta*[ 1+2*epsilon^2 -(1+epsilon); -(1+epsilon) 2+epsilon^2];

%% Compute the New Sum
y0 = X_obs;
Sum_e_NEW = (norm(Rbar * (xhat_Givens - Xbar)))^2 + (X_obs(1) - H(1,:)*xhat_Givens)^2 + (X_obs(2) - H(2,:)*xhat_Givens)^2;

%% Return
P_Givens    = P2_Givens;
P_True      = P2_True;
X_Givens      = xhat_Givens;
