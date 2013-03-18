function BI = triad(B1, B2, I1, I2)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
%                           triad.m
% Author:   Zach Dischner
% Date:     Feb 28, 2013
% 
% Usage:
%   BI = triad(b1, b2, i1, i2)
%
% Description:  Given two body frame observations, and their Intertial
%               representations, determine the rotation matrix [BI] 
%               between the two
% 
% Inputs:  b1 => body vector 1       (most accurate)
%          b2 => body vector 2       (least accurate)
%          I1 => inertial vector 1   (most accurate)
%          I2 => inertial vector 2   (least accurate)
%
% Outputs: BI => Standard DCM representing the rotation between body
%                and intertial frames
% 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Preallocate
IT = zeros(3,3);    % DCM between triad and inertial frame
BT = zeros(3,3);    % DCM between triad and body frame


%% Triad baby!
% First triad frame axis is colinear with the 'better' observation
BT(:,1) = B1; IT(:,1) = I1; 

% Second triad frame axis is the norm'd cross product between observed axes
BT(:,2) = cross(B1,B2)/norm(cross(B1,B2));
IT(:,2) = cross(I1,I2)/norm(cross(I1,I2));

% Third triad frame axes is simply the cross product of the two 
BT(:,3) = cross(BT(:,1),BT(:,2));
IT(:,3) = cross(IT(:,1),IT(:,2));

% Calculate rotation matrix:
BI = BT*IT';

end