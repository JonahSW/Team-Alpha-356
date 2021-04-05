function [fid, oev1, jdpp] = readoe2(filename)

% read orbital elements data file

% required by encke.m and cowell2.m

% input

%  filename = name of orbital element data file

% output

%  fid    = file id
%  jdpp   = julian date of perihelion passage
%  oev(1) = semimajor axis
%  oev(2) = orbital eccentricity
%  oev(3) = orbital inclination
%  oev(4) = argument of perigee
%  oev(5) = right ascension of the ascending node
%  oev(6) = true anomaly

% NOTE: all angular elements are returned in radians

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dtr = pi / 180;

% open data file

fid = fopen(filename, 'r');

% check for file open error

if (fid == -1)
    clc; home;
    fprintf('\n\n error: cannot find this file!!');
    keycheck;
    return;
end

% read 27 lines of data file

for i = 1:1:27
    cline = fgetl(fid);
    switch i
        case 3
            tl = size(cline);

            ci = findstr(cline, ',');

            % extract month, day and year

            month = str2double(cline(1:ci(1)-1));

            day = str2double(cline(ci(1)+1:ci(2)-1));

            year = str2double(cline(ci(2)+1:tl(2)));
        case 7
            tl = size(cline);

            ci = findstr(cline, ',');

            % extract hours, minutes and seconds

            uthr = str2double(cline(1:ci(1)-1));

            utmin = str2double(cline(ci(1)+1:ci(2)-1));

            utsec = str2double(cline(ci(2)+1:tl(2)));

            jdpp = julian(month, day, year) ...
                + uthr / 24 + utmin / 1440 + utsec / 86400;
        case 11
            oev1(1) = str2double(cline);
        case 15
            oev1(2) = str2double(cline);
        case 19
            oev1(3) = dtr * str2double(cline);
        case 23
            oev1(4) = dtr * str2double(cline);
        case 27
            oev1(5) = dtr * str2double(cline);
    end
end

status = fclose(fid);

