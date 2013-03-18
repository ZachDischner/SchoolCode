function e2 = subEuler231(e,e1)

% subEuler231(E,E1)
%
%	E2 = subEuler231(E,E1) computes the relative
%	(2-3-1) Euler angle vector from E1 to E.
%

C = Euler2312C(e);
C1 = Euler2312C(e1);
C2 = C*C1';
e2 = C2Euler231(C2);
