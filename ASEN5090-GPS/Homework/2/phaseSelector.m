function [p1,p2] = phaseSelector(PRNNumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Sept 16 2013
%      +-------------------------+------------------------\    __      ____
%  {=  |     Zach Dischner,      | Email:   zach.dischner  \___| \____/ ___\___|
%  {=  |  Aerospace Engineering  |         @colorado.edu   /   ===---======)--->
%      +-------------------------+------------------------/                    |
%
%     Function  =>
%
%     Purpose   => Find the correct "phases" or indices you need to
%                   generate a PRN code for say a GPS satellite
%                  Made for ASEN 5090
%
%     Input     => PRNNumber : Number of the satellite you are searching
%                              for the indices for
%
%     Output    => p1,p2     : Two indices for using in a PRN generator
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

phaseTable=[   2, 6;
            3 , 7;
            4 , 8;
            5 , 9;
            1 , 9;
            2 , 10;
            1 , 8;
            2 , 9;
            3 , 10;
            2 , 3;
            3 , 4;
            5 , 6;
            6 , 7;
            7 , 8;
            8 , 9;
            9 , 10;
            1 , 4;
            2 , 5;
            3 , 6;
            4 , 7;
            5 , 8;
            6 , 9;
            1 , 3;
            4 , 6;
            5 , 7;
            6 , 8;
            7 , 9;
            8 , 10;
            1 , 6;
            2 , 7;
            3 , 8;
            4 , 9;
            5 , 10;
            4 , 10;
            1 , 7;
            2 , 8;
            4 , 10];
p1 = phaseTable(PRNNumber,1);
p2 = phaseTable(PRNNumber,2);
