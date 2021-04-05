%jjs280
%04/4/2021
%This code is used to estimate the deltaV required for a low thrust trajectory between circular coplanar orbits

%% Constants
%Sun
mu_sun = 1.32712440018e20; %[m^3/s^2]
AU = 149597870.7e3;% [m]
%Earth
mu_earth = 3.986004418e14; %[m^3/s^2]
earth_ephemeris = read_ephemeris(2);
a_earth = mean(earth_ephemeris(3,:))*AU;
in_earth = mean(earth_ephemeris(4,:));
r_earth = 6378000;% [m]
r_earth_SOI = 0.929e9;% [m]
earth_axial_tilt = 23.43928*pi/180;% [rad]
%Ceres
mu_ceres = 6.26325e10; %[m^3/s^2]
r_ceres = 469730;% [m]
r_ceres_SOI = 7709455.15;% [m]
ceres_ephemeris = read_ephemeris(4);
a_ceres = mean(ceres_ephemeris(3,:))*AU;
in_ceres = mean(ceres_ephemeris(4,:));

%Alpha Spacecraft
LEO_assembly_altitude = 400e3;% [m]
LEO_capture_altitude = 600e3;% [m]
LCO_altitude = 500e3;%[m]
launch_inclination = 28.573*(pi/180);

%DeltaV for leaving Earth's SOI (1)
r0 = LEO_assembly_altitude + r_earth;
r = r_earth_SOI;
deltaV1 = spiral(mu_earth,r0,r);

%DeltaV for traveling to Ceres (2)
r0 = a_earth;
r = a_ceres;
deltaV2 = spiral(mu_sun,r0,r);

%DeltaV for capturing in Ceres SOI (3)
r0 = r_ceres_SOI;
r = LCO_altitude + r_ceres;
deltaV3 = spiral(mu_ceres,r0,r);

%DeltaV for leaving Ceres' SOI (4)
r0 = LCO_altitude + r_ceres;
r = r_ceres_SOI;
deltaV4 = spiral(mu_ceres,r0,r);

%DeltaV for traveling to Earth (5)
r0 = a_ceres;
r = a_earth;
deltaV5 = spiral(mu_sun,r0,r);

%DeltaV for capturing in Earth SOI (6)
r0 = r_earth_SOI;
r = LEO_capture_altitude + r_earth;
deltaV6 = spiral(mu_earth,r0,r);

%DeltaV for changing inclination in Earth Orbit


%DeltaV for changing inclination in Heliocentric Orbit


%% Results
deltaV_total = (deltaV1+deltaV2+deltaV3+deltaV4+deltaV5+deltaV6)*1e-3;
disp(['The total mission deltaV is ',num2str(deltaV_total),' km/s.']);

%%
function deltaV = spiral(mu,r0,r)
    deltaV = abs(sqrt(mu/r0) - sqrt(mu/r));
end

function deltaV = inclination()
    
end