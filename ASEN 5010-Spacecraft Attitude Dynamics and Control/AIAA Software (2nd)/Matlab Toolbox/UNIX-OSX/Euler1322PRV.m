function	q = Euler1322PRV(e)

% Euler1322PRV(E)
%
%	Q = Euler1322PRV(E) translates the (1-3-2) Euler
%	angle vector E into the principal rotation vector Q.
%

q = EP2PRV(Euler1322EP(e));
