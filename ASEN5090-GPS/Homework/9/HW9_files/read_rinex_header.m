%**************************************************************************
% function [ fid, rec_xyz, observables ] = read_rinex_header( file_name )
% 
% This file takes a RINEX header file and parses it for
% the approximate location of the receiver in XYZ
% coordinates (meters) and for the number and types of 
% observables contained in the corresponding RINEX file.
% 
% Arguments:
% 
% file_name - the name of the header file.  Note, you must put the 
%             file name in '' marks when calling header.m.
%             Example:  read_header('header.txt')
% 
% Returns:
% 
% fid - the file ID number used for the header file.
% rec_xyz - the approximate receiver location XYZ in meters.
% observabes - the observables contained in the file, this is a
%              1 x n cell where n = number of observables.
%              
% By Jan Weiss
% Last Updated 2-17-03
%**************************************************************************

function [ fid, rec_xyz, observables ] = read_rinex_header( file_name )

% Initialize vars
observables = {};
rec_xyz = [ NaN NaN NaN ];

% Assign a file ID and open the given header file.
fid=fopen(file_name);

% If the file does not exist, scream bloody murder!
if fid == -1
    display('Error!  Header file does not exist.');
else
    
    % Set up a flag for when the header file is done.
    end_of_header=0;
    
    % Get the first line of the file.
    current_line = fgetl(fid);
    
    % If not at the end of the file, search for the desired information.
    while end_of_header ~= 1
        
        % Search for the approximate receiver location line.
        if strfind(current_line,'APPROX POSITION XYZ')
            
            % Read xyz coordinates into a matrix.
            rec_xyz = sscanf(current_line,'%f');
        end
        
        % Search for the number/types of observables line.
        if strfind(current_line,'# / TYPES OF OBSERV')
            
            % Read the non-white space characters into a temp variable.
            [observables_temp] = sscanf(current_line,'%s');            
            
            % Read the number of observables space and then create
            % a matrix containing the actual observables.
            for ii = 1:str2num(observables_temp(1))                
                observables{ii} = observables_temp( 2*ii : 2*ii+1 );
            end
          
        end
        
                  
        % Get the next line of the header file.
        current_line = fgetl(fid);
        
        %Check if this line is at the end of the header file.
        if strfind(current_line,'END OF HEADER')
            end_of_header=1;
        end
        
    end
end