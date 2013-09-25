function [az,el,range] = ASEN5090_ecef2azelrange(r_sat,r_site,latgd,lon)
%==========================================================================
% [az,el,range] = ecef2azelrange(r_sat,r_site,latgd,lon)
%
% Calculates the azimuth, elevation, and range of a satellite with respect
%  to an observation site.
%
% Author: Ben K. Bradley
% Date: 11/15/2010
% Homework changes made by Zach Dischner and Nick Truesdale, 9/24/13
%
% INPUT:         Description                                         Units
%  r_sat      - position of satellite in ECEF frame                 [x y z]
%  r_site     - position of observing site in ECEF frame            [x y z]
%  latgd      - geodetic latitude of observation site          [-90,90] deg
%  lon        - longitude of observation site     [-180,180] or [0,360] deg
%
% OUTPUT:
%  az         - azimuth (degrees clockwise from North)          [0,360] deg
%  el         - elevation (degrees up from horizon)            [-90,90] deg
%  range      - distance from observation site to satellite
%==========================================================================


% Range vector, ECEF
pECF = r_sat - r_site;

% Convert to SEZ
pSEZ = rotmat(2, pi/2 - latgd*pi/180)*rotmat(3, lon*180/pi)*pECF';

% Range
range = norm(pSEZ);

% Azimuth
az = wrapTo360( atan2d(pSEZ(2), -pSEZ(1)) );

% Elevation
el = asind(pSEZ(3)/range);


end %function

function B = rotmat(dim, theta, A)
%
% File: rotation_matrix.m
% Author: Nick Truesdale, my partner on this homework assignment. This is a
% helpful rotation matrix function for use. 
% Date: 9/21/2012
%
% Description: This function calculates a rotation matrix in either two or
% three dimensions. DIM = 0 yields a 2D matrix, while DIM = 1,2 or 3 yields
% one of the three 3D matrices. THETA is the angle of rotation. If A is
% given as a vector or array, B will be the rotated vector/array. If A is
% unspecified, B will be the rotation matrix.
%
% Note: The size of A must match the given dimension. If DIM = 0, A must be
% a 2x2 or 2x1 matrix. If DIM = 1,2,3, A must be a 3x3 or 3x1 matrix.

% Parse inputs
mat_dim = logical(dim) + 2;
if(nargin < 3), A = []; end

% Sine and cosine
S = sin(theta);
C = cos(theta);

% Calculate matrix
switch dim
    case 0
        B = [C, S;
            -S, C];
    case 1
        B = [1  0 0;
            0  C S;
            0 -S C];
    case 2
        B = [C 0 -S;
            0 1  0;
            S 0  C];
    case 3
        B = [C  S  0;
            -S  C  0;
            0  0  1];
    otherwise
        error('Dim must be 0,1,2 or 3, denoting a 2D, X, Y or Z rotation')
end

% If A is supplied
if(~isempty(A))
    if(isequal(size(A), [mat_dim, mat_dim]) || isequal(size(A), [mat_dim, 1]))
        B = B*A;
    else
        error(['A must be a column vector or square matrix with dimension ',...
            '2 for DIM = 0, or dimension 3 for DIM = 1,2,3'])
    end
end
end