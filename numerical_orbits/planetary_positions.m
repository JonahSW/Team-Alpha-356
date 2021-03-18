%jjs280
%03/09/2021
%Generates heliocentric position vectors using planetary ephemeris data from JPL HORIZONS Tool
%Ephemeris data given per day
%Linearly interpolate between days for smaller timesteps

%Plots orbits if plot == 1
function plotions = planetary_positions(plotme, animateme)

%Target Body = Earth [Geocenter]
r_earth = position_vector(2);

%Target Body = Ceres
r_ceres = position_vector(4);

%Target Body = Mars
r_mars = position_vector(3);

%Target Body = Venus
r_venus = position_vector(1);

%Target Body = Jupiter
r_jupiter = position_vector(5);

plotions = vertcat(r_venus,r_earth,r_mars,r_ceres,r_jupiter);

% Plot Results
if animateme == 1
   animate(r_ceres); 
end

if plotme == 1
    figure()
    polarplot(r_venus(1,:),r_venus(2,:), r_earth(1,:),r_earth(2,:), r_mars(1,:),r_mars(2,:), r_ceres(1,:),r_ceres(2,:), r_jupiter(1,:),r_jupiter(2,:));
    
    %Generate cartesian coords
    [venus_x, venus_y, venus_z] = sphere2cart(r_venus(2,:),r_venus(1,:),r_venus(3,:));
    [earth_x, earth_y, earth_z] = sphere2cart(r_earth(2,:),r_earth(1,:),r_earth(3,:));
    [mars_x, mars_y, mars_z] = sphere2cart(r_mars(2,:),r_mars(1,:),r_mars(3,:));
    [ceres_x, ceres_y, ceres_z] = sphere2cart(r_ceres(2,:),r_ceres(1,:),r_ceres(3,:));
    [jupiter_x, jupiter_y, jupiter_z] = sphere2cart(r_jupiter(2,:),r_jupiter(1,:),r_jupiter(3,:));
    
    figure()
    plot3(0,0,0,'*y');%Plot the sun
    hold on
    grid minor
    axis equal
    plot3(venus_x, venus_y, venus_z);
    plot3(earth_x, earth_y, earth_z);
    plot3(mars_x, mars_y, mars_z);
    plot3(ceres_x, ceres_y, ceres_z);
    plot3(jupiter_x, jupiter_y, jupiter_z);
    
end

end

%% Supporting Functions

%Function for calculating the position vector of a planet
%Position relative to the heliocentric coordinate system in spherical coordinates
function r = position_vector(body)
    %Generate position vector from ephemeris data
    matrix = read_ephemeris(body);
    eccentricity = matrix(1,:)';
    true_anomaly = matrix(2,:)';
    semi_major_axis = matrix(3,:)';
    inclination = matrix(4,:)';
    
    len = length(matrix(1,:));
    
    r = zeros(3,len);
    for i=1:1:len
        a = semi_major_axis(i);
        e = eccentricity(i);
        theta = true_anomaly(i);% Convert degrees to radians
        phi = inclination(i)*sin(theta);% Use inclination to find position vector, convert degrees to radians
        r(1,i) = theta;
        r(2,i) = (a*(1-e^2)) / (1+e*cos(theta));% Calculate radius from true anomaly and semi-major axis
        r(3,i) = phi;
    end
    
end

%Converts speherical to cartesian coordinates for plotting in 3d
function [x, y, z] = sphere2cart(r,theta,phi)
    x = r .* cos(phi) .* cos(theta);
    y = r .* cos(phi) .* sin(theta);
    z = r .* sin(phi);
end

%Function for animating orbits
%In Progress
function animate(position_vector)
    figure()
    polarplot(position_vector(1,1),position_vector(2,1));
    hold on
    for i = 2:1:length(position_vector)
       polarplot(position_vector(1,i),position_vector(2,i),'o');       
    end
    hold off
end

%Horizons Tool:
%https://ssd.jpl.nasa.gov/horizons.cgi?CGISESSID=83291ab9cbf5c73d25e0aef2b6067c20#results
%Inputs:
%Ephemeris Type = Elements
%Center = Sun (body center)
%Time Span: Start=2040-01-01, Stop=2050-01-01, Step=1 d