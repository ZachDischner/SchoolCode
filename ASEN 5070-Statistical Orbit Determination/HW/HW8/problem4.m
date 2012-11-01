function problem4()


%% a - Derive P1
problem4()
syms(e)
H       = [1 2*e;1 3*e];
R       = eye(2,2);
P0bar   = (1/e^2)*eye(2,2);
P1      = @(e) inv( transpose(H)*inv(R) * H + P0bar);

