function q = OLAE(Vb, Vi, w)

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           OLAE.m
% Author:   Zach Dischner
% Date:     April 20, 2013
% 
% Usage:
%   q = OLAE( Vb, Vi ) (W is optional)
%
% Description:  Uses the Optimal Linear Attitude Estimation (OLAE)
%               algorithm to find a least squares solution for the DCM 
%               between the B and N frames. 
% 
%               Lamens ==> Solves Wabba's problem in a super cool linear 
%                          way. Finds [BN] 
% 
% Inputs:  Vb  ==> Body frame vectors for different measurements
%                       Looks like [ Vb1 Vb2 ...Vbn ]
%                       where Vbn is a column vector
%          Vi  ==> Corresponding inertial frame vectors for different measurements
%                       Looks like [ Vi1 Vi2 ...Vin ]
%                       where Vin is a column vector
% %        W   ==> Optional weights to apply to each measurement
%                   [W1 W2...] becomes | W1 0 ... |
%                                      |  0 W2 ...|
% Outputs: q   ==> Classical Rodriguez param representing rotation between
%                   B and N frames
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Assemble Weighting Matrix
numvecs = size(Vb,2);
if nargin < 3
    w = ones(numvecs,1);
end

W = zeros(3*length(w), 3*length(w/norm(w)));
for ii = 1:length(w)
    W(1+3*(ii-1):3+3*(ii-1),1+3*(ii-1):3+3*(ii-1)) = w(ii)*eye(3,3);
end


%% Assemble 'd' and 'S' matrices
d = Vb - Vi;
d = reshape(d,numel(d),1);

s = Vb + Vi;
S = zeros(3*numvecs,3);

for ii=1:numvecs
    S(1+3*(ii-1):3+3*(ii-1),1:3) = tilde(s(:,ii));
end

%% Solve it brah!
q = (S'*W*S)\(S'*W*d);


    




