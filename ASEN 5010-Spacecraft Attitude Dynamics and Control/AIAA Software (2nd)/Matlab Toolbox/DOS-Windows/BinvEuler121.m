function B = BinvEuler121(q)

% BinvEuler121(Q)
%
%	B = BinvEuler121(Q) returns the 3x3 matrix which relates 
%	the derivative of the (1-2-1) Euler angle vector Q to the 
%	body angular velocity vector w.  
%	
%		w = [B(Q)]^(-1) dQ/dt
%	


s2 = sin(q(2));
c2 = cos(q(2));
s3 = sin(q(3));
c3 = cos(q(3));

B(1,1) = c2;
B(1,2) = 0;
B(1,3) = 1;
B(2,1) = s2*s3;
B(2,2) = c3;
B(2,3) = 0;
B(3,1) = s2*c3;
B(3,2) = -s3;
B(3,3) = 0;
