function w=intRigidEOM(w_init,t,I)

%% Rotation about each PA

w = zeros(length(t),length(w_init));
w(1,:)= w_init;
for ii=2:length(t)
   % Find the Derivative
   wprime = RigidEOMDot(t(ii),w(ii-1,:),I);
   
   % General linear integration:
   %    x_(n+1) = x_(n) + x'*delta_t
   w(ii,:) = w(ii-1,:) + wprime*(t(ii)-t(ii-1));
   
end




