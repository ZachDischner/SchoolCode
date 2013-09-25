function PRN = PRNGenerator(PRNNumber,epoch)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Sept 16 2013
%      +-------------------------+------------------------\    __      ____
%  {=  |     Zach Dischner,      | Email:   zach.dischner  \___| \____/ ___\___|
%  {=  |  Aerospace Engineering  |         @colorado.edu   /   ===---======)--->
%      +-------------------------+------------------------/                    |
%
%     Function  => PRNGenerator
%
%     Purpose   => Generate PRN Code for a given satellite number
%
%     Input     => PRNNumber : Number of the satellite you are searching
%                              for the indices for
%
%     Output    => PRN     : 1024 bit PRN code (as an array)
%
%     Procedure =>
%
%                 __...____________________          ,
%                `(\ [ ===NCC-1700===--|__|) ___..--"_`--.._____
%                  `"""""""""""""""""| |""` [_""_-___________"_/
%                                    | |   /..../`'-._.-'`
%                                ____| |__/::..'_
%                               |\ ".`"` '_____//\
%                               `"'-.   """""  \\/
%                                    `""""""""""`
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==1
    epoch=0;
end
%% Allocate Initial Variables
% Two shift registers
n   = 10;
len = 2^n-1;
G1  = ones(1,n);
G2  = G1;
G1Mask  = [0 0 1 0 0 0 0 0 0 1];
G2Mask  = [0 1 1 0 0 1 0 1 1 1];

% Output code
PRN = zeros(1,len);

%% Get 'Phase' Indices
[G2idx(1),G2idx(2)] = phaseSelector(PRNNumber);
G1idx = [3,10];

%% Create PRN Code
for ii=1:len+epoch
   G2out = xor( G2(G2idx(1)), G2(G2idx(2))) ;
   if ii > epoch
        PRN(ii-epoch) = xor(G1(end),G2out);
   end
   G1 = [ mod(sum(G1.*G1Mask),2), G1(1:n-1) ];
   G2 = [ mod(sum(G2.*G2Mask),2), G2(1:n-1) ];
end



