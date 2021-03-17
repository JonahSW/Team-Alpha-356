%jjs280
%03/16/2021
%Implements a function that returns the ephemeris data for a given input
%1 = venus, 2 = earth, 3 = mars, 4 = ceres, 5 = jupiter

%Output:
% 1->eccentricity
% 2->true anomaly
% 3->semi-major axis
% 4->inclination

function data = read_ephemeris(a)
    if a == 1
        matrix = readmatrix('ephemeris_data\venus_ephemeris.csv');
        inclination = matrix(:,3)*(pi/180);
        eccentricity = matrix(:,1);
        true_anomaly = matrix(:,9)*(pi/180);
        semi_major_axis = matrix(:,10);
        data = vertcat(eccentricity',true_anomaly',semi_major_axis',inclination');
    elseif a == 2
        matrix = readmatrix('ephemeris_data\earth_ephemeris.csv');
        inclination = matrix(:,3)*(pi/180);
        eccentricity = matrix(:,1);
        true_anomaly = matrix(:,9)*(pi/180);
        semi_major_axis = matrix(:,10);
        data = vertcat(eccentricity',true_anomaly',semi_major_axis',inclination');
    elseif a == 3
        matrix = readmatrix('ephemeris_data\mars_ephemeris.csv');
        inclination = matrix(:,3)*(pi/180);
        eccentricity = matrix(:,1);
        true_anomaly = matrix(:,9)*(pi/180);
        semi_major_axis = matrix(:,10);
        data = vertcat(eccentricity',true_anomaly',semi_major_axis',inclination');
    elseif a == 4
        matrix = readmatrix('ephemeris_data\ceres_ephemeris.csv');
        inclination = matrix(:,3)*(pi/180);
        eccentricity = matrix(:,1);
        true_anomaly = matrix(:,9)*(pi/180);
        semi_major_axis = matrix(:,10);
        data = vertcat(eccentricity',true_anomaly',semi_major_axis',inclination');
    elseif a == 5
        matrix = readmatrix('ephemeris_data\jupiter_ephemeris.csv');
        inclination = matrix(:,3)*(pi/180);
        eccentricity = matrix(:,1);
        true_anomaly = matrix(:,9)*(pi/180);
        semi_major_axis = matrix(:,10);
        data = vertcat(eccentricity',true_anomaly',semi_major_axis',inclination');
    else
        disp('Could not locate ephemeris data.');
        disp('');
        disp('Please enter 1, 2, 3, 4, or 5 to select venus, earth, mars, ceres, or jupiter');
    end
end