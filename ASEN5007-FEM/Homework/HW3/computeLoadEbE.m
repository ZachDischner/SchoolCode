function [f1] = computeLoadEbE(q0, qi, qj, L1, L2)
% Loads acting on f1 based on the surrounding elements
% 
%                  qj      
%       qi         |         |
%  q0    |         |         
%  |_L1__|___L2____|
%        f1       f2   
f1 = L1/6*(q0 + 2*qi) + L2/6*(2*qi + qj);
