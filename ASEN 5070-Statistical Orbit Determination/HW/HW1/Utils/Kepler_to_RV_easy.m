function[Rxyz, Vxyz] = Kepler_to_RV_easy(a, e, i, Omega, w, nu, mu)

p = a*(1 - e^2);                    % semiparameter (km)

%%%  Position Coordinates in Perifocal Coordinate System
x  = (p*cos(nu)) / (1 + e*cos(nu)); % x-coordinate (km)
y  = (p*sin(nu)) / (1 + e*cos(nu)); % y-coordinate (km)
z  = 0;                             % z-coordinate (km)
vx = -(mu/p)^(1/2) * sin(nu);       % velocity in x (km/s)
vy = (mu/p)^(1/2) * (e + cos(nu));  % velocity in y (km/s)
vz = 0;                             % velocity in z (km/s)

%%%  Transformation Matrix (3 Rotations)  %%%
ROT = [cos(Omega)*cos(w)-sin(Omega)*sin(w)*cos(i) ...
                (-1)*cos(Omega)*sin(w)-sin(Omega)*cos(w)*cos(i) ...
                            sin(Omega)*sin(i); ...
       sin(Omega)*cos(w)+cos(Omega)*sin(w)*cos(i) ...
                (-1)*sin(Omega)*sin(w)+cos(Omega)*cos(w)*cos(i) ...
                            (-1)*cos(Omega)*sin(i); ...
       sin(w)*sin(i)  cos(w)*sin(i)  cos(i)];

%%%  Transforming Perifocal -> xyz  %%%
Rxyz = ROT*[x y z]';
Vxyz = ROT*[vx vy vz]';