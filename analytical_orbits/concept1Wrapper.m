%jjs280
%03/17/2021
%Calculates the delta V and mass ratios required for a concept 1 mission
%Assumes circular orbits and applies patched conics
close all
clc
clear

%Vehicle Properties
dry_mass = 90e3;
payload_mass = 3e3;
tank_mass = 9e3;
num_tanks = 7;
Isp = 950;

%Inputs:
ceres_orbit_altitude = 70e3;
earth_orbit_altitude = 400e3;

%Solar System Constants
%Earth
earth_ephemeris = read_ephemeris(2);
a_earth = mean(earth_ephemeris(3,:));
in_earth = mean(earth_ephemeris(4,:));
earth_tilt = 23.43928*(pi/180);
launch_inclination = 28.573*(pi/180);
earth_radius = 6378000;
mu_earth = 3.986004418e14;% Earth Gravitational Parameter [m^3/sec^2]
    
%Ceres
ceres_ephemeris = read_ephemeris(4);
a_ceres = mean(ceres_ephemeris(3,:));
in_ceres = mean(ceres_ephemeris(4,:));
ceres_tilt = 4*(pi/180);
ceres_radius = 469730;
mu_ceres = 6.26325e10;% Ceres Gravitational Parameter [m^3/sec^2]

%Earth Departure (plane change + departure burn)
delta = (earth_tilt-launch_inclination)*(pi/180)+0.9*(in_earth+in_ceres);% 90% of heliocentric plane change is done at Earth
a_start = earth_radius+earth_orbit_altitude;
[Mratio_1, deltaV_1, psi_1] = hohmann_departure_burn(a_earth, a_ceres, mu_earth, delta, a_start, Isp);

%Ceres Capture (plane change + capture burn)
delta = 0.1*(in_earth+in_ceres);% 10% of heliocentric plane change is done at Ceres
a_end = ceres_radius+ceres_orbit_altitude;
[Mratio_2, deltaV_2, psi_2] = hohmann_capture_burn(a_earth, a_ceres, mu_ceres, delta, a_end, Isp);

%Ceres Departure (plane change + departure burn)
delta = 0.1*(in_earth+in_ceres);% 10% of heliocentric plane change is done at Ceres
a_start = ceres_radius+ceres_orbit_altitude;
[Mratio_3, deltaV_3, psi_3] = hohmann_departure_burn(a_ceres, a_earth, mu_ceres, delta, a_start, Isp);

%Earth Capture (capture burn)
delta = 0.9*(in_earth+in_ceres);% 90% of heliocentric plane change is done at Earth
a_end = earth_radius+earth_orbit_altitude;
[Mratio_4, deltaV_4, psi_4] = hohmann_capture_burn(a_earth, a_ceres, mu_earth, delta, a_end, Isp);%Switching a_earth and a_ceres

%Spacecraft Mass
%After burn4
m4 = (dry_mass + 1*tank_mass);
%After burn3
m3 = (m4 + 2*tank_mass)*Mratio_4;
%After burn2
m2 = (m3 + 2*tank_mass + payload_mass)*Mratio_3;
%After burn1
m1 = (m2 + 3*tank_mass)*Mratio_2;
%Initial
m0 = (m1 + 4*tank_mass)*Mratio_1;

%Mission Duration
duration_outbound = hohmann_duration(a_earth, a_ceres);
wait = wait_time(a_earth, a_ceres);
duration_inbound = hohmann_duration(a_earth, a_ceres);

%% Plot Results
figure()
plot(1:1:4,[deltaV_1,deltaV_2,deltaV_3,deltaV_4],'o');
title('Delta V For Each Burn');
xlabel('Burn Number');
ylabel('Delta V (m/s)');

figure()
plot(1:1:5,[m0,m1,m2,m3,m4],'o');
title('Mass');
ylabel('Mass (kg)');

disp(['The outbound trajectory duration is ',num2str(duration_outbound),' days'])
disp(['The wait time at Ceres is ',num2str(wait),' days'])
disp(['The inbound trajectory duration is ',num2str(duration_outbound),' days'])
disp(['The total duration is ',num2str(duration_inbound + wait + duration_outbound),' days'])