function plotBinaryWave(X,t)

% Simple function, plots binary wave (X) in a pretty form
% optional t input for x range
if nargin==1
    t = linspace(1,length(X),length(X));
end
stairs(t,X); 
axis([t(1),t(end),-1,2]);
