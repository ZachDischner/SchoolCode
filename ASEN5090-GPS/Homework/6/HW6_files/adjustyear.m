function year = adjust_year(yy)
% P. Axelrad - function to correct 2 digit year to the right value

yr=1900*(yy<100); % Convert 00-99 to 1900-1999
yr = yr + 100*(yy<30); % Convert 00-30 to 2000-2030

year = yy+yr; 
