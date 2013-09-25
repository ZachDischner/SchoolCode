%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        quickFEM.m
% Author      : Zach Dischner
% Date        : 9/24/2013
% Description : quick Matlab checker script for Problem 4, HW 3
%
%                 __...____________________          ,
%                `(\ [ ===NCC-1700===--|__|) ___..--"_`--.._____
%                  `"""""""""""""""""| |""` [_""_-___________"_/
%                                    | |   /..../`'-._.-'`
%                                ____| |__/::..'_
%                               |\ ".`"` '_____//\
%                               `"'-.   """""  \\/
%                                    `""""""""""`
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%% Part 1, the X forces
d = 180-[180 136 84 36 0];
pp=62.4;
p = pp*d;
p = [0 p 0];
d = [0 d 0];
L = (d - [d(end) d(1:end-1)]); L(L<0) = 0;
for ii = 2:length(p)-1
    Fx(ii-1) = computeLoadEbE(p(ii-1),p(ii), p(ii+1),L(ii),L(ii+1));
end
fprintf('Computed X Node loads are:\n')
Fx

%% Part 2, the Y forces
d = [0 70 210 350 490];
p = -1*ones(length(d),1)'*180*pp;
p = [0 p 0];
d = [0 d 0];
L  = d - [d(end) d(1:end-1)]; L(L<0) = 0;
for ii = 2:length(p)-1
    Fy(ii-1) = computeLoadEbE(p(ii-1),p(ii), p(ii+1),L(ii),L(ii+1));
end
fprintf('Computed Y Node loads are:\n')
Fy



