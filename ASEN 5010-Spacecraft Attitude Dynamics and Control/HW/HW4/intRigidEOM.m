function X=intRigidEOM(X_init,t,I,returnPrime)

%% Rotation about each PA

% X is a column vector with [w ; 321angles]

X = zeros(length(t),length(X_init));
Xprime = X;
X(1,:)= X_init;

for ii=2:length(t)
   % Find the Derivative
   Xprime(ii,:) = RigidEOMDot(t(ii),X(ii-1,:),I);
   
   % General linear integration:
   %    x_(n+1) = x_(n) + x'*delta_t
   X(ii,:) = X(ii-1,:) + Xprime(ii,:)*(t(ii)-t(ii-1));
   
   
end

if exist('returnPrime','var')
    X=Xprime;
end




