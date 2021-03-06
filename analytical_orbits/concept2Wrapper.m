%jjs280
%03/17/2021
%Calculates the delta V and mass ratios required for a concept 2 mission
%Assumes circular orbits and applies patched conics
close all
clc
clear

%Vehicle Properties
dry_mass_crewed = 68.5e3;
dry_mass_tug = 178.58e3;
payload_mass_crewed = 500;
payload_mass_tug = 4e3;%Payload that can be abandoned
tank_mass_NTP = 8.55e3;
tank_mass_NEP = 2.43e3;
num_NTP_tanks_crewed = 2;
num_NTP_tanks_tug = 2;
num_NEP_tanks_tug = 2;
boil_off_mass_transfer = 6e3;%32.06e3;%Boil off for 1/2 Hohmann transfer
boil_off_mass_LEO = 7e3;%55e3;
boil_off_mass_Ceres = 1e3;%39e3;

Isp_NTP = 950;
Isp_NEP = 4000;
thrust_NEP = 80.1;

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

%Spacecraft Mass
%Crewed Vessel
%After burn4
m5 = (dry_mass_crewed + 1*tank_mass_NTP + payload_mass_crewed);
%After burn3
m4 = (m5 + 1*tank_mass_NTP + boil_off_mass_transfer)*Mratio_4;
%After Resupply
m3 = (m4 + boil_off_mass_Ceres)*Mratio_3;
%After burn2
m2 = (dry_mass_crewed + 1*tank_mass_NTP + payload_mass_crewed + boil_off_mass_transfer);
%After burn1
m1 = (m2 + 1*tank_mass_NTP)*Mratio_2;
%Initial
m0 = (m1 + boil_off_mass_LEO)*Mratio_1;

%Fuel Tug
Mratio_NEP_1 = exp((deltaV_1+deltaV_2)/(g0*Isp_NEP));
%End
m2_NEP = (dry_mass_tug + 2*tank_mass_NEP);
%After burn 1
m1_NEP = (dry_mass_tug + payload_mass_tug + 2*tank_mass_NTP + boil_off_mass_Ceres + ((m4-m5) + (m3-m4)));
%Initial
m0_NEP = (m1_NEP + boil_off_mass_transfer)*Mratio_NEP_1;

%Mission Duration
%Hohmann Transfer
duration_outbound = hohmann_duration(a_earth, a_ceres);
wait = wait_time(a_earth, a_ceres);
duration_inbound = hohmann_duration(a_earth, a_ceres);
%Low Thrust Estimation - Outbound
duration_NEP = low_thrust_estimator((deltaV_1+deltaV_2), Isp_NEP, m0_NEP, thrust_NEP);

%% Plot Results
figure()
plot(1:1:4,[Mratio_1,Mratio_2,Mratio_3,Mratio_4]);
title('Mass Ratio');

figure()
plot(1:1:6,[m0,m1,m2,m3,m4,m5],'o');
hold on
grid minor
title('Mass');
plot(1:1:6,[m0_NEP,m1_NEP,m2_NEP,m2_NEP,m2_NEP,m2_NEP],'o');
xline(1,'-','LEO Assembly Mass');
xline(2,'-','After Earth Departure Burn');
xline(3,'-','After Ceres Capture Burn');
xline(4,'-','After Fuel Transfer');
xline(5,'-','After Ceres Departure Burn');
xline(6,'-','After Earth Capture Burn','LabelHorizontalAlignment','left');
ylabel('Mass (kg)');
legend('Crewed Vehicle','Fuel Tug');

disp(['The crewed vessel outbound trajectory duration is ',num2str(duration_outbound),' days'])
disp(['The crewed vessel wait time at Ceres is ',num2str(wait),' days'])
disp(['The crewed vessel inbound trajectory duration is ',num2str(duration_inbound),' days'])
disp(['The crewed vessel total mission duration is ',num2str(duration_inbound + wait + duration_outbound),' days, or ',num2str((duration_inbound + wait + duration_outbound)/365.25),' years'])

disp(['The fuel tug outbound trajectory duration is ',num2str(duration_NEP),' days, or ',num2str((duration_NEP)/365.25),' years'])