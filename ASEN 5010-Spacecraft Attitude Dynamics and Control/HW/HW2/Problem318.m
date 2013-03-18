clc;clear all;close all
fprintf('Problem 3.18:\n\n\n')
%% Create Symbolic Variables
syms a0 a1 a2 a3 b0 b1 b2 b3

% a# is q'
% b# is q''

%% Build DCMs and Multiple
C1 = EP2C([a0 a1 a2 a3]);

C2 = EP2C([b0 b1 b2 b3]);

C = C2*C1;

%% Use Sheperd's equations

% Tell it what 'one' is!!!
% Can look at the trace term and see that this is a common element to
% simplify

one = (a0^2 + a1^2 + a2^2 + a3^2) * (b0^2 + b1^2 + b2^2 + b3^2); 
% Good trick, thank you kevin dinkel


B0 = simplify(sqrt(1/4*(one + trace(C))));

B1 = simplify(sqrt(1/4*(one+2*C(1,1)-trace(C))));

B2 = simplify(sqrt(1/4*(one+2*C(2,2)-trace(C))));

B3 = simplify(sqrt(1/4*(one+2*C(3,3)-trace(C))));

%% Output
fprintf('\nB0 = ');pretty(B0)
fprintf('\nB1 = ');pretty(B1)
fprintf('\nB2 = ');pretty(B2)
fprintf('\nB3 = ');pretty(B3)

fprintf('\n\n Note, there is sign ambiguity on a symbolic square root')
fprintf('\n So Matlab doesnt know how to treat it. That is why there')
fprintf('\n a [...^2)^(1/2)] term in the final answers. They are irrelevant')
fprintf('\n\n Also note that B0 is negative of what is in the book')
fprintf('\n Since quaternions are a dual set, this is just the opposite')
fprintf('\n set for a rotation\n\n')