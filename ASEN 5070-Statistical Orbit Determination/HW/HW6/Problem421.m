function [a b c] = Problem421()


%% Givens
Pbar0   = [4 0 0;0 2 0;0 0 1];
Xbar0   = [1 1 1]';
Yt1     = 2;                    
xt1     = 2;
varEy   = 1;        % Variance of E(Y(t1)). That is, Sigma^2(eY)
R       = varEy;


%% (a) - Compute Xhat(t0) using the batch processor
Phi = @(t) [1 t t^2/2;0 1 t;0 0 1];

 
G = [1 0 0];
Phi1 = Phi(1);

Htilde1 = [1 0 0]; %dG/dx   maps X to Y. Y=(X1)=2

H1 = Htilde1*Phi1;

Y1 = 2;

Lamda = inv(Pbar0) + H1'*inv(R)*H1;
N     = (Pbar0)\Xbar0 + H1'*inv(R)*Y1;

Xhat0 = inv(Lamda)*N;

a=Xhat0;

format rat
%% (b) - Sequential
% use many of the same variables as before. Essentially use Kalman stuff
% with H, not Htilde1 to map back to zero. 


K0    = Pbar0*H1'*inv(H1*Pbar0*H1' + R);

Xhat0 = Xbar0 + K0*(Y1 - H1*Xbar0);  % Estimating X0 with original guess, and a statistical measure of 

b     = Xhat0;


%% (c) - Phi
% New Covariance Matrix
Pbar1 = Phi1*Pbar0*Phi1';
% New Xbar
Xbar1 = Phi1 * Xbar0;


% Kalman Gain
K1    = Pbar1*Htilde1'*inv(Htilde1*Pbar1*Htilde1' + R);

Xhat1 = Xbar1 + K1*(Y1 - Htilde1*Xbar1);



c = Phi1\Xhat1;





end