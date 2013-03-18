function s = MRPswitch(q,s2)

% MRPswitch
%
%	S = MRPswitch(Q,s2) checks to see if norm(Q) is larger than s2.
%	If yes, then the MRP vector Q is mapped to its shadow set.
%
%   % PASS IN ROWS DUDE!!!!!

q2 = dot(q,q);   % Smarter than q'q
if (q2>s2*s2)
	s = -q/q2;
else
	s = q;
end
