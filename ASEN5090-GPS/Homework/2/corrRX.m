function Rx= corrRX(x,n)
% x is the function
% n is the shift
len=length(x);
Rx = 1/len * sum(x.*[x(end-n:end),x(1:end-n-1)]);
