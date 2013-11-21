%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                        getIonFreePseudorange.m
% Author      : Zach Dischner
% Date        : 10/24/2013
% Description : Convert a date type object into [GPS_Weeks, GPS_SOW] time
%
%
%                 __...____________________          ,
%                `(\ [ ===NCC-1700===--|__|) ___..--"_`--.._____
%                  `"""""""""""""""""| |""` [_""_-___________"_/
%                                    | |   /..../`'-._.-'`
%                                ____| |__/::..'_
%                               |\ ".`"` '_____//\
%                               `"'-.   """""  \\/
%                                    `""""""""""`
% Inputs        : obs_data - observed dataset from rined
%                 observables - list of strings representing the observable
%                   values in the observation file
%
% Outputs       : rho_3 - ion free pseudorange
%                       * Uses P2 and P1 if it can, otherwise it used P2
%                       and C1
%
% TODOS         : If I wan't to use this, I'm going to need to pass in a
%                 lot more variables. Maybe a structure of a bunch of stuff
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function rho_3 = getIonFreePseudorange(obs_data, observables) 

P2 = obs_data.data(data_idx,P2_col);

if sum(strcmp(observables,'P1'))>0
    %------Retrieve P1 as 2nd freq pseudorange
    PorC1 = obs_data.data(data_idx,P1_col);
else
    %------Retrieve C1 as 2nd freq pseudorange
    PorC1 = obs_data.data(data_idx,C1_col);
end

%------Use ONLY ionosphere free. Uncomment linesto use what data we have
if P2 == 0 || PorC1 == 0
    rho_3 = -1;
    fprintf('CANT USE ION-FREE PSEUDORANGE!!!!, returning -1');
else
    rho_3 = 2.5457*PorC1-1.5457*P2;
end
