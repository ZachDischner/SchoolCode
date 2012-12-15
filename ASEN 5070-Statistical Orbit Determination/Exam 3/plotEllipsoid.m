% Plot an ellipsoid given an orthonormal, right handed
% transformation matrix, R and the semi - axis, semi
%
% For the Stat. O.D. project R is made up of the eigenvectors
% of the upper 3x3 portion of the covariance matrix.  semi
% contains sigma_x, sigma_y, sigma_z in a column vector.  

function plotEllipsoid(R,semi)

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