%jjs280
%03/17/2021
%Calculates the delta V and mass ratios required for a concept 2 mission
%Assumes circular orbits and applies patched conics
close all
clc

%Vehicle Properties
dry_mass_crewed = 90e3;
dry_mass_tug = 60e3;
payload_mass_crewed = 1e3;
payload_mass_tug = 2e3;
tank_mass_NTP = 9e3;
tank_mass_NEP = 4.5e3;
num_NTP_tanks_crewed = 3;
num_NTP_tanks_tug = 3;
num_NEP_tanks_tug = 2;

Isp_NTP = 950;
Isp_NEP = 4000;
thrust_NEP = 20;

%Inputs:
ceres_orbit_altitude = 70e3;
earth_orbit_altitude = 400e3;

%Solar System Constants
%Earth
g0 = 9.807;
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
[Mratio_1, deltaV_1, psi_1] = hohmann_departure_burn(a_earth, a_ceres, mu_earth, delta, a_start, Isp_NTP);
%Ceres Capture (plane change + capture burn)
delta = 0.1*(in_earth+in_ceres);% 10% of heliocentric plane change is done at Ceres
a_end = ceres_radius+ceres_orbit_altitude;
[Mratio_2, deltaV_2, psi_2] = hohmann_capture_burn(a_earth, a_ceres, mu_ceres, delta, a_end, Isp_NTP);
%Ceres Departure (plane change + departure burn)
delta = 0.1*(in_earth+in_ceres);% 10% of heliocentric plane change is done at Ceres
a_start = ceres_radius+ceres_orbit_altitude;
[Mratio_3, deltaV_3, psi_3] = hohmann_departure_burn(a_ceres, a_earth, mu_ceres, delta, a_start, Isp_NTP);
%Earth Capture (capture burn)
delta = 0.9*(in_earth+in_ceres);% 90% of heliocentric plane change is done at Earth
a_end = earth_radius+earth_orbit_altitude;
[Mratio_4, deltaV_4, psi_4] = hohmann_capture_burn(a_earth, a_ceres, mu_earth, delta, a_end, Isp_NTP);

%Low Thrust Estimation (Outbound fuel tug)
duration_NEP = low_thrust_estimator((deltaV_1+deltaV_2), Isp_NEP, dry_mass_tug, thrust_NEP);
Mratio_NEP_1 = exp((deltaV_1+deltaV_2)/(g0*Isp_NEP));

%Spacecraft Mass
%Crewed Vessel
%After burn4
m5 = (dry_mass_crewed + 1*tank_mass_NTP + payload_mass_crewed);
%After burn3
m4 = (m5 + 2*tank_mass_NTP)*Mratio_4;
%After Resupply
m3 = (m4)*Mratio_3;
%After burn2
m2 = (dry_mass_crewed + 1*tank_mass_NTP + payload_mass_crewed);
%After burn1
m1 = (m2 + 2*tank_mass_NTP)*Mratio_2;
%Initial
m0 = (m1)*Mratio_1;

%Fuel Tug
%End
m2_NEP = (dry_mass_tug + 2*tank_mass_NEP + payload_mass_tug);
%After burn 1
m1_NEP = (dry_mass_tug + 2*tank_mass_NEP + payload_mass_tug + 3*tank_mass_NTP);
%Initial
m0_NEP = (m1)*Mratio_NEP_1;

%Mission Duration
%Hohmann Transfer
duration_outbound = hohmann_duration(a_earth, a_ceres);
wait = wait_time(a_earth, a_ceres);
duration_inbound = hohmann_duration(a_earth, a_ceres);
%Low Thrust Estimation - Outbound
duration_NEP = low_thrust_estimator((deltaV_1+deltaV_2), Isp, dry_mass, thrust_NEP);

%Plot Results
figure()
plot(1:1:4,[Mratio_1,Mratio_2,Mratio_3,Mratio_4]);
title('Mass Ratio');

figure()
plot(1:1:6,[m0,m1,m2,m3,m4,m5],'o');
hold on
title('Mass');
plot(1:1:6,[m0_NEP,m1_NEP,m2_NEP,m2_NEP,m2_NEP,m2_NEP],'o');

disp(['The crewed vessel outbound trajectory duration is ',num2str(duration_outbound),' days'])
disp(['The crewed vessel wait time at Ceres is ',num2str(wait),' days'])
disp(['The crewed vessel inbound trajectory duration is ',num2str(duration_inbound),' days'])
disp(['The crewed vessel total mission duration is ',num2str(duration_inbound + wait + duration_outbound),' days'])

disp(['The fuel tug outbound trajectory duration is ',num2str(duration_NEP),' days'])