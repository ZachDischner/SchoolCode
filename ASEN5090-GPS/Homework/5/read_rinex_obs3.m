function [ rinexv3 ] = read_rinex_obs(fname, PRN_list, maxlines)
%*******************************************************
% function [ rinexv3 ] = read_rinex(fname, PRN_list)
%
% DESCRIPTION:
%
%     This function reads a RINEX format GPS data
%     file and returns the data in an array.
%
% ARGUMENTS:
%
%     fname (str) - RINEX file
%     PRN_list (opt) - vector with PRNs to process, useful
%                      for ignoring some PRN data.
%
% OUTPUT:
%
%     rinexv3 - rinex data (v3 struct)
%
% CALLED BY:
%
%
% FUNCTIONS CALLED:
%
%     read_rinex_header.m
%
% MODIFICATIONS:
%
%     XX-XX-03  :  Jan Weiss - Original
%     03-14-06  :  Jan Weiss - Cleanup
%               :  See SVN log for further modifications
%     10-26-09  :  P. Axelrad to get rid of conversion of phase to m and
%     validity checking for ASEN 5090
%     10-14-11  :  P. Axelrad - put in check for GPS satellites only
%     10-17-11  :  P. Axelrad - corrected time conversion call
%               :  put in check for comment lines and added C2 data type
%     11-2-12  :  P. Axelrad - now correctly handles files with more than 12
%     satellites tracked and Glonass, and speeds up reading of large files.
%
% Colorado Center for Astrodynamics Research
% Copyright 2006 University of Colorado, Boulder
%*******************************************************
blocksize = 100000;

% Initialize variables
rinex_data = [];
line_count = 1;


% Read header

[ fid, rec_xyz, observables ] = read_rinex_header(fname);
num_obs = length(observables);

% Reserve array
if num_obs > 5
    r_data = zeros(blocksize,13);
else
    r_data = zeros(blocksize,8);
end
icount = 1;


% Status
disp([ 'Parsing RINEX file ' fname ]);
if nargin >= 2 & PRN_list
    disp([ 'Returning data for PRNs: ' num2str(PRN_list) ]);
else
    disp('Returning data for all PRNs');
end

if nargin < 3
    maxlines = 4000000;
end

% Get the first line of the observations.
current_line = fgetl(fid);

% If not at the end of the file, search for the desired information.
while (current_line ~= -1) & (line_count < maxlines)
    
    
    %Check if this line is a comment line.
    if isempty(strfind(current_line,'COMMENT') )
        
        
        % Check event flag
        event_flag = str2num(current_line(27:29));
        if event_flag == 0
            yr = adjustyear(str2num(current_line(2:3)));
            
            % Get the time for this data epoch.
            current_time = [ yr ; str2num(current_line(5:6)) ; ...
                str2num(current_line(8:9)) ; str2num(current_line(11:12)) ; ...
                str2num(current_line(14:15)) ; str2num(current_line(16:26)) ]';
            
            [~, gpssec, gpswk] = gpsvec2gpstow(current_time);
            
            
            
            % How many SV's are there?
            current_num_sv = str2num(current_line(30:32));
            current_prn=zeros(current_num_sv,1);
            sat_type=repmat('X',current_num_sv,1);
            
            % Figure out which PRN's there are.
            
            num_sat_lines = ceil(current_num_sv/12)-1; % read extra line(s) if too many satellite for one line
            for i = 1:num_sat_lines
                temp = fgetl(fid);
                current_line = [current_line temp(33:end)];
            end
            for ii=1:current_num_sv
                sat_type(ii) = current_line(30 + 3*ii);
                current_prn(ii) = str2num(current_line(31 + 3*ii : 32 + 3*ii));
            end
            
            
            % Get the data for all SV's in this epoch.
            for ii=1:current_num_sv
                
                % Get the next line.
                current_line = fgetl(fid);
                line_count = line_count + 1;
                if rem(line_count, 1000) == 0
                    disp([ 'Read ', num2str(line_count) ' lines' ]);
                end
                
                
                % Check the length of the line and pad it with zeros to
                % make sure it is 80 characters long.
                current_line = check_rinex_line_length(current_line);
                
                % Get the observables on this line.
                current_obs = [ str2num(current_line(1:14)) ; str2num(current_line(17:30)) ; ...
                    str2num(current_line(33:46)) ; str2num(current_line(49:62)) ; str2num(current_line(65:78)) ];
                
                % If there are > 5 observables, read another line to get the rest of the observables for this SV.
                if num_obs > 5
                    
                    % Get the next line.
                    current_line = fgetl(fid);
                    line_count = line_count + 1;
                    
                    
                    % Check the length of the line and pad it with zeros to
                    % make sure it is 80 characters long.
                    current_line = check_rinex_line_length(current_line);
                    
                    % Append the data in this line to the data from previous line.
                    current_obs = [ current_obs ; str2num(current_line(1:14)) ; ...
                        str2num(current_line(17:30)) ; str2num(current_line(33:46)) ; ...
                        str2num(current_line(49:62)) ; str2num(current_line(65:78)) ];
                    
                end  % if num_obs > 5
                
                
                
                
                % Glonass add 100 to PRN, Galileo add 200 to PRN, Compass add 300
                % to PRN
                switch sat_type(ii)
                    case 'R'
                        current_prn(ii)=current_prn(ii)+100;
                    case 'E'
                        current_prn(ii)=current_prn(ii)+200;
                    case 'C'
                        current_prn(ii)=current_prn(ii)+300;
                    case ~('G' | ' ')
                        current_prn(ii)=0;
                        continue
                end
                
                % Format the data for this PRN as Date/Time, PRN, Observations.
                current_data = [ gpswk, gpssec, current_prn(ii) , current_obs'];
                
                % Keep only data for the specified PRNs
                if nargin >= 2 & PRN_list & isempty(find(PRN_list == current_prn(ii)))
                    continue
                end
                
                
                %Append to the master rinex data file.
                %rinex_data = [ rinex_data ; current_data ];
                r_data(icount,:)= current_data;
                if icount == blocksize
                    rinex_data = [rinex_data; r_data];
                    r_data = zeros(size(r_data));
                    icount = 1;
                else
                    icount = icount+1;
                end
                
                
                
                
            end  % for ii=1:current_num_sv
            
        end % for event flag
    end % for comment line
    % Get the next line.
    current_line = fgetl(fid);
    line_count = line_count + 1;
    
    
end  % while current_line ~= -1
rinex_data = [rinex_data; r_data];
i = find(rinex_data(:,2) == 0);
rinex_data(i,:) =[];

size(rinex_data)
rinexv3.data = rinex_data;
clear rinex_data

% Define columns
rinexv3 = define_cols(rinexv3, observables);

% Get rid of zero data entries introduced by
% line padding when more than 5 obs are present
% rinexv3.data = rinexv3.data(:,1:3+num_obs);


% Status
disp([ 'Total lines: ', num2str(line_count) ]);
disp('Finished.');
disp(' ');



% =========================================================================
function rinex = define_cols(rinex, observables)

% Defaults
rinex.col.WEEK = 1;
rinex.col.TOW = 2;
rinex.col.PRN = 3;

col_offset = 3;

for ii=1:length(observables)
    
    switch observables{ii}
        case {'L1'}
            rinex.col.L1 = ii + col_offset;
        case {'L2'}
            rinex.col.L2 = ii + col_offset;
        case {'LA'}
            rinex.col.LA = ii + col_offset;
        case {'P1'}
            rinex.col.P1 = ii + col_offset;
        case {'P2'}
            rinex.col.P2 = ii + col_offset;
        case {'C1'}
            rinex.col.C1 = ii + col_offset;
        case {'C2'}
            rinex.col.C2 = ii + col_offset;
        case {'S1'}
            rinex.col.S1 = ii + col_offset;
        case {'S2'}
            rinex.col.S2 = ii + col_offset;
        case {'SA'}
            rinex.col.SA = ii + col_offset;
    end  % switch
    
end  % for ii=1:length(observables)


function [ current_line ] = check_rinex_line_length(current_line)

if length(current_line) < 80
    
    add_spaces = 80 - length(current_line);
    
    for j = 1 : add_spaces
        
        current_line = [ current_line , '0' ];
        
    end
    
end

% Check if there are any blanks in the data and put a zero there.
current_line = strrep(current_line,' ', '0');




