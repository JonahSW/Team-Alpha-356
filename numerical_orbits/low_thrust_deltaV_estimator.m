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
LEO_assembly_altitude = 120e3;% [m]
LEO_capture_altitude = 1000e3;% [m]
LCO_altitude = 700e3;%[m]
launch_inclination = 28.573*(pi/180);
Isp = 8510;% [s]
g0 = 9.807;% [m/s^2]
m_dry = 220e3;% [kg]
m_dot = 1.343e-3;%[kg/s]

%DeltaV for leaving Earth's SOI (1)
r0 = LEO_assembly_altitude + r_earth;
r = r_earth_SOI;
v1 = sqrt(mu_earth/r0);
v2 = sqrt(mu_earth/r);
delta_i = (launch_inclination - earth_axial_tilt)*180/pi;% plane change for leaving Earth, deg
deltaV1 = spiral(mu_earth,r0,r) + inclination_change(v1,v2,delta_i);

%DeltaV for traveling to Ceres (2)
r0 = a_earth;
r = a_ceres;
v1 = sqrt(mu_sun/r0);
v2 = sqrt(mu_sun/r);
delta_i = (in_earth + in_ceres)*180/pi;% plane around sun, deg
deltaV2 = spiral(mu_sun,r0,r) + inclination_change(v1,v2,delta_i);

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
v1 = sqrt(mu_sun/r0);
v2 = sqrt(mu_sun/r);
delta_i = (in_earth + in_ceres)*180/pi;% plane around sun, deg
deltaV5 = spiral(mu_sun,r0,r) + inclination_change(v1,v2,delta_i);

%DeltaV for capturing in Earth SOI (6)
r0 = r_earth_SOI;
r = LEO_capture_altitude + r_earth;
deltaV6 = spiral(mu_earth,r0,r);

%% Results
deltaV_total = (deltaV1+deltaV2+deltaV3+deltaV4+deltaV5+deltaV6);
disp(['The total mission deltaV is ',num2str(deltaV_total*1e-3),' km/s.']);

%Spacecraft Mass
mratio = exp(deltaV_total/(Isp*g0));
m_wet = m_dry*mratio;
m_propellant = m_wet-m_dry;
duration = m_propellant/m_dot;
disp(['The total mission mass ratio is ',num2str(mratio),'']);
disp(['The total spacecraft wet mass is ',num2str(m_wet),' kg.']);
disp(['The total propellant mass is ',num2str(m_propellant),' kg.']);
disp(['The total thrusting duration is ',num2str(duration/(3600*24)),' days.']);

%%
function deltaV = spiral(mu,r0,r)
    deltaV = abs(sqrt(mu/r0) - sqrt(mu/r));
end

function deltaV = inclination_change(v1,v2,delta_i)
    deltaV = sqrt(v1^2 + v2^2 - 2*v1*v2*cos((delta_i*pi/180)*pi/2));
end