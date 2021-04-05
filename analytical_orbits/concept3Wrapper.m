%jjs280
%03/17/2021
%Calculates the delta V and mass ratios required for a concept 3 mission
%Assumes circular orbits and applies patched conics
close all
clc
clear

%Vehicle Properties - UPDATE dry_mass and thrust from mass budget calculations and iterate until close
dry_mass = 161113.8;%Includes tanks, minus payload and kick stages
payload_mass = 5470.9;
Ceres_sample_mass = 10;
kick_stage_mass = 2500;
boil_off_mass_transfer = 1e2;%11.02e3;%Boil off for 1/2 Hohmann transfer
boil_off_mass_LEO = 1e2;%8.56e3;
boil_off_mass_Ceres = 1e2;%7.59e3;
Isp = 4000;
thrust = 63.2;%N
%Inputs:
ceres_orbit_altitude = 200e3;
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

%% Estimate deltaV costs using hohmann transfers
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
[Mratio_4, deltaV_4, psi_4] = hohmann_capture_burn(a_earth, a_ceres, mu_earth, delta, a_end, Isp);

%Estimate delta V from low thrust trajectory - rough guesses, unused
duration = 2*365.25*24*3600;
wet_mass = 300e3;
deltaV_outbound = thrust*duration/((dry_mass+wet_mass)*2/3);
deltaV_inbound = thrust*duration/((dry_mass+wet_mass)*1/3);

%Low Thrust Estimation - Outbound
Mratio_NEP_1 = exp((deltaV_1+deltaV_2)/(g0*Isp));%Mass ratio for outbound trajectory
%Low Thrust Estimation - Inbound
Mratio_NEP_2 = exp((deltaV_3+deltaV_4)/(g0*Isp));%Mass ratio for inbound trajectory

%Spacecraft Mass
%After burn 2
m3 = (dry_mass + kick_stage_mass);
%At Ceres Departure
m2 = (m3 + boil_off_mass_transfer + kick_stage_mass + Ceres_sample_mass)*Mratio_NEP_2;
%After burn 1
m1 = (m2 + payload_mass + boil_off_mass_Ceres + kick_stage_mass);
%Before burn 1
m0 = (m1)*Mratio_NEP_1 + boil_off_mass_transfer + boil_off_mass_LEO + kick_stage_mass;

%Propellant mass
m_propellant = (m2-m3)+(m0-m1);

%Wait Time:
duration1 = low_thrust_estimator((deltaV_1+deltaV_2), Isp, m0, thrust);
duration2 = low_thrust_estimator((deltaV_3+deltaV_4), Isp, m2, thrust);
wait = 180;%Currently a very rough estimate

%% Plot Results
figure()
plot(1:1:4,[m0,m1,m2,m3],'o');
grid minor
xline(1,'-','LEO Assembly Mass');
xline(2,'-','After Outbound Low Thrust Trajectory');
xline(3,'-','At Ceres Departure');
xline(4,'-','After Inbound Low Thrust Trajectory','LabelHorizontalAlignment','left');
title('Mass');
ylabel('Mass (kg)');

disp(['The outbound trajectory duration is ',num2str(duration1),' days'])
disp(['The wait time at Ceres is ',num2str(wait),' days'])
disp(['The inbound trajectory duration is ',num2str(duration2),' days'])
disp(['The total duration is ',num2str(duration1 + wait + duration2),' days, or ',num2str((duration1 + wait + duration2)/365.25),' years'])
disp(['The spacecraft wet mass is ',num2str(m0),' kg'])
disp(['The required propellant mass is ',num2str(m_propellant),' kg'])