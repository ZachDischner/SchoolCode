% Plot an ellipsoid given an orthonormal, right handed
% transformation matrix, R and the semi - axis, semi
%
% For the Stat. O.D. project R is made up of the eigenvectors
% of the upper 3x3 portion of the covariance matrix.  semi
% contains sigma_x, sigma_y, sigma_z in a column vector.  

function plotEllipsoidP(P)

[evecs,evals] = eig(P);
semi(1) = sqrt(evals(1,1));
semi(2) = sqrt(evals(2,2));
semi(3) = sqrt(evals(3,3));
semi = sort(semi);
semi = semi([3,2,1]);

[x,y,z] = sphere(20);

x = x * semi(1); 
y = y * semi(2);
z = z * semi(3);

[mm,nn] = size(x);

C = (R * [x(:) y(:) z(:)]')';

x = reshape(C(:,1),mm,nn);
y = reshape(C(:,2),mm,nn);
z = reshape(C(:,3),mm,nn);

surf(x,y,z)
shading flat
axis equal
xlabel('X');ylabel('Y');zlabel('Z')