function [a b c] = Problem420()

%% Givens
y           = [1 2 1]';
Ee          = 0;
R           = [1/2 0 0;0 1 0;0 0 1];   % E(ee')
H           = [1 1 1]';
xbar        = 2;
var_xbar    = 2;                       % sigma^2(xbar)

%% Part A: Find xhat with batch processing algorithm
W   = inv(R);
Wbar=var_xbar;
xhat = (H'*W*H + Wbar)\(H'*W*y+Wbar*xbar);

a    = xhat;



%% Part B: Find Standard Deviation of the Estimation Error with xhat
epsilon = y-H*xhat;
b       = std(epsilon);


%% Part C: Find Estimation Error
c       = epsilon;



end