%jjs280
%03/09/2021
%Generates heliocentric position vectors using planetary ephemeris data from JPL HORIZONS Tool
%Ephemeris data given per day
%Linearly interpolate between days for smaller timesteps

%Plots orbits if plot == 1
function ephemerides = ephemerides()

%Target Body = Earth [Geocenter]
inputFile = 'earth_ephemeris.csv';
r_earth = position_vector(inputFile);

%Target Body = Ceres
inputFile = 'ceres_ephemeris.csv';
r_ceres = position_vector(inputFile);

%Target Body = Mars
inputFile = 'mars_ephemeris.csv';
r_mars = position_vector(inputFile);

%Target Body = Venus
inputFile = 'venus_ephemeris.csv';
r_venus = position_vector(inputFile);

%Target Body = Jupiter
inputFile = 'jupiter_ephemeris.csv';
r_jupiter = position_vector(inputFile);

ephemerides = vertcat(r_venus,r_earth,r_mars,r_ceres,r_jupiter);

% Plot Results
plot = 1;
if plot == 1
    figure(1)
    polarplot(r_venus(1,:),r_venus(2,:), r_earth(1,:),r_earth(2,:), r_mars(1,:),r_mars(2,:), r_ceres(1,:),r_ceres(2,:), r_jupiter(1,:),r_jupiter(2,:));

    figure (2)
    [cart_venus_x, cart_venus_y, cart_venus_z] = sph2cart(r_venus(2,:)',r_venus(3,:)',r_venus(1,:)');
    plot3(cart_venus_x, cart_venus_y, cart_venus_z)
end

end

%Supporting Function
function r = position_vector(inputFile)
    %Generate position vector from ephemeris data
    matrix = readmatrix(inputFile);
    inclination = matrix(:,3);
    eccentricity = matrix(:,1);
    true_anomaly = matrix(:,9);
    semi_major_axis = matrix(:,10);
    
    len = length(matrix(:,1));
    
    r = zeros(2,len);
    for i=1:1:len
        a = semi_major_axis(i);
        e = eccentricity(i);
        theta = true_anomaly(i)*(pi/180);% Convert degrees to radians
        phi = inclination(i)*(pi/180);% Convert degrees to radians
        r(1,i) = theta;
        r(2,i) = (a*(1-e^2)) / (1+e*cos(theta));% Calculate radius from true anomaly
        r(3,i) = phi;
    end
    
end

%Horizons Tool:
%https://ssd.jpl.nasa.gov/horizons.cgi?CGISESSID=83291ab9cbf5c73d25e0aef2b6067c20#results
%Inputs:
%Ephemeris Type = Elements
%Center = Sun (body center)
%Time Span: Start=2040-01-01, Stop=2050-01-01, Step=1 d