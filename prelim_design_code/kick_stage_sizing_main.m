%jjs280
%03/31/2021
%This code is used to calculate the performance of the chemical kick stage
%Calculate deltaV to get to Earth SOI and required mass ratio with a chemical kick stage
clc;
%clear;
close all;


%Hohmann Transfer from assembly orbit to intermediate orbit:
mu_earth = 3.986e14; %[m^3/s^2]
R_earth = 6.371e6; %[m]
initial_altitude = 8e5;% [m]
R_im = (R_earth + initial_altitude):1e6:2e8; %[m]
a1 = R_earth + initial_altitude;
a2 = R_im;
g0 = 9.807; %[m/s]
%Calculate deltaV
deltaV_1 = sqrt(mu_earth./a1).*(sqrt((2.*a2)./(a1+a2))-1);
deltaV_2 = sqrt(mu_earth./a2).*(1-sqrt((2.*a1)./(a1+a2)));
%Kick Stage Sizing (Biprop, Full Hohmann to MEO)
Isp = 400; %[s]
spacecraft_dry_mass = 300e3; %[kg]
Mratio = exp((deltaV_1+deltaV_2)/(Isp*g0));
propellant_mass = spacecraft_dry_mass*(Mratio - 1);
%disp(['Propellant mass for a Hohmann Transfer from LEO to MEO is: ',num2str(propellant_mass),' kg']);
figure()
plot(a2.*1e-3,deltaV_1,a2.*1e-3,deltaV_2,a2.*1e-3,(deltaV_2+deltaV_1));
title('DeltaV vs. Final Orbit Radius (m)');
legend('burn 1', 'burn 2', 'total');
figure()
plot(a2*1e-3,Mratio);
title('Mass Ratio vs. Final Orbit Altitude');

%Boost Stage Sizing (Big Boost)
%120000kg Biprop, 110000 kg propellant, Isp = 400s;
Isp = 400;% [s]
spacecraft_dry_mass = 3.5e5+1e4; %[kg]
solid_propellant_mass = 1e5;% [kg]
deltaV_kick = Isp*g0*log((spacecraft_dry_mass+solid_propellant_mass)/spacecraft_dry_mass);
disp(['Delta V for the boost stage is: ',num2str(deltaV_kick),' m/s']);

%Kick Stage Sizing (Partial kick)
%4x 6000kg SRs, 5000 kg propellant each, Isp = 300s, burn 4 at once;
Isp = 300;% [s]
spacecraft_dry_mass = 3.5528e5+2.4e4; %[kg]
solid_propellant_mass = 4*5e3;% [kg]
deltaV_kick = Isp*g0*log((spacecraft_dry_mass+solid_propellant_mass)/spacecraft_dry_mass);
disp(['Delta V for the kick stage is: ',num2str(deltaV_kick),' m/s']);


%Hohmann Transfer from assembly orbit to Earth SOI:
%{
%CONCLUSION: Chemical kick stage cannot be used to carry s/c to Earth SOI, it would need >400 mT of propellant
mu_earth = 3.986e14; %[m^3/s^2]
R_earth = 6.371e6; %[m]
R_earth_SOI = 0.929e9; %[m]
initial_altitude = 600e3;% [m]
a1 = R_earth + initial_altitude;
a2 = R_earth_SOI;
g0 = 9.807; %[m/s]
%Calculate deltaV
deltaV_1 = sqrt(mu_earth/a1)*(sqrt((2*a2)/(a1+a2))-1);
deltaV_2 = sqrt(mu_earth/a2)*(1-sqrt((2*a1)/(a1+a2)));
%Kick Stage Sizing (Biprop, Full Hohmann to SOI)
Isp = 400; %[s]
spacecraft_dry_mass = 300e3; %[kg]
Mratio = exp((deltaV_1+deltaV_2)/(Isp*g0));
propellant_mass = spacecraft_dry_mass*(Mratio - 1);
disp(['Propellant mass for a Hohmann Transfer from LEO to Earth SOI is: ',num2str(propellant_mass),' kg']);
%}