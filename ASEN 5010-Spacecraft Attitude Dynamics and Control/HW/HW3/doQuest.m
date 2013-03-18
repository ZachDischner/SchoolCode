function q = doQuest(vB,vI,weights)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           doQuest.m
% Author:   Zach Dischner
% Date:     Feb 28, 2013
% 
% Usage:
%   BI = doQuest(vB,vI)
%
% Description:  Computes the rotation between two frames using a least
%               squares fit between frame observations
% 
% Inputs:  vB => Array of body observations [vb1  vb2 ...]
%          vI => Array of Inertial observations [vi1  vi2 ;...]
%          weights => (optional) weighting matrix for each set of
%          observations
%
% Outputs: Q => Euler parameters representing the rotation DCM [BI]
%               Order is: [q1i,q2j,q3k,q4]
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Setup
if ~exist('weights','var') 
    weights = ones(1,size(vB,2));     % uniform weighting
end

% Initial guess for lambda
lam = sum(weights);
% f = @(s,k) det(K-s*eye(4,4));
%% Get the parameters of K, not K itself though

    B = bsxfun(@times,weights,vB)*vI';
    S = B + B';
    Z = [  B(2,3)-B(3,2);
           B(3,1)-B(1,3);
           B(1,2)-B(2,1)];

    sig = trace(B);
    K = [S-sig*eye(3),Z;Z',sig];
    

p = inv((lam +sig)*eye(3,3)-S)*Z;

%% Convert CRP to Quaternions

q = 1/sqrt(1+p'*p)*[p;1];