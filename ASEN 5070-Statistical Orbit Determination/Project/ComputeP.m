function [P2]=ComputeP(H,P1bar,R,Xstar0,Xobs,method)
% Compute P2 (or the trace therof) for various different algorithms. 

% Input is error value, and 'tr', indicator to return the trace of each 
%       calculation



% Looks like    y=Hx+e
H1      = H(1,:);
H2      = H(2,:);

if method == 'kalman'
%% b - P2 Using Conventional Kalman Filter
K1 = P1bar*H1'*inv(H1*P1bar*H1' + R(1,1));
P1 = (eye(size(K1*H1)) - K1*H1)*P1bar;

% Now again for P2
K2          = P1*H2'*inv(H2*P1*H2' + R(2,2));
P2   = (eye(size(K2*H2)) - K2*H2)*P1;
%------------------------------------------------------------

end

if method == 'joseph'
%% c - P2 Using Joseph Formulation
% Branching off of the Kalman
K1          = P1bar*H1'*inv(H1*P1bar*H1' + R(1,1));
P1          = (eye(size(K1*H1)) - K1*H1)*P1bar*(eye(size(K1*H1))-K1*H1)' + K1*R(1,1)*K1';
K2          = P1*H2'*inv(H2*P1*H2' + R(2,2));

P2   = (eye(size(K2*H2)) - K2*H2)*P1*(eye(size(K2*H2))-K2*H2)' + K2*R(2,2)*K2';
%------------------------------------------------------------
end

if method == 'potter'
%% d - P2 Using Plotter Algorithm (5.7.17) 
[R,p]=chol(P1bar);
if p ==0
   W1bar   = (chol(P1bar))'; 
else
    W1bar   = (chol(P1bar,'lower'));
end

Ftilde  = W1bar'*H1';
alpha   = inv(Ftilde'*Ftilde + R(1,1));
gamma   = 1./(1 + sqrt(R(1,1)*alpha));
K       = alpha*W1bar*Ftilde;
% xhat    = xbar + K*([v1;v2]);
W1      = W1bar - gamma*K*Ftilde';
P1      = W1*W1';

% Now for obs 2
W2bar   = chol(P1)';
Ftilde  = W2bar'*H2';
alpha   = inv(Ftilde'*Ftilde + R(2,2));
gamma   = 1./(1 + sqrt(R(2,2)*alpha));
K       = alpha*W2bar*Ftilde;
% xhat    = xbar + K*([v1;v2]);
W2      = W2bar - gamma*K*Ftilde';

P2=  W2*W2';

end


if method == 'givens'
    h = H;
    P0bar = P1bar;
    xbar=Xstar0';
    y=Xobs;
    % Givens Algorithm
%------------------------------------------------------
% Find a-priori stuff for  Givens
% Choleski decompisition

[R,p]=chol(P0bar);
if p ==0
   Sbar   = (chol(P0bar))'; 
else
    Sbar   = (chol(P0bar,'lower'));
end

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


m       = length(Xobs);
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


P2 = inv(U)*inv(D)*inv(U)';
end



          