function Xprime = diffMRPLin(t,X,I,P,K)

Xprime = [   zeros(3,3)    ,    eye(3,3)/4     ; ...
            -inv(I)*K  ,   -inv(I)*P]*X';
