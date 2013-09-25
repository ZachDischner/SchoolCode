function [lag, Rxy]= cyc_corr2(x,y)
%function [lag, Rxy] = cyc_corr2(x,y)
% inputs
%   y is the reference signal
%   x is the one you want to compare to the reference signal
% outputs
%  lag (for the x-axis)
%  Rxy = correlations (for the y-axis)
%
% kristine's comments: I got this code from Penny, and have made a few small changes.

% original code
Y = [y y y]';
n = length(x);
for i = 1:2*n % sweep the lag from -?n to + n
    Rxy(i) = x*Y(i:i+n-1);
    lag(i) = i-1-n;
end

% Let's normalize it
Rxy = Rxy/n;
% take off the first point because that one just is a repeat
Rxy = Rxy(2:end);
lag = lag(2:end);