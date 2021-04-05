% sdf33
% updated 4.3.21

%% Ceres Craft Maneuvers
clear; clc; close all;

% selection initial altitude
% assume all spacecrafts begin at LCO on elliptic plane
% must be calculated for each craft
% finds delta V, mass, mass ratios, duration for each transfer

%% Missions
% Polar Study Missions: Polar Lander 1 & Polar Lander 2
% General Surface Study Mission: General Lander 1 & 2, Seismometer Lander 1 & 2
% Surface Mapping Study: Mapping Orbiter
% Communications Array Mission: Comms Orbiter 1, 2, 3, & 4

%% Initial Parameters
% parameters
mu = 62630000000; % m3/s2, Ceres gravitation parameter
axial = 0.06981317008; % rad, axial tilt
radius = 469730; % m, average radius of surface
g = 9.807; % m/s2, acceleration of Earth
soi = 7.7e7; % m

% low Ceres orbit
alt = 100000; % altitude of LCO
a2 = radius + alt; % m, radius of LCO



%% Polar Study Mission
fprintf('Polar Study Mission: Polar Lander 1 & 2\n\n')

% landing site = north pole
angdeg = 90+axial;
ang = deg2rad(angdeg);

% low Ceres orbit --> surface (hohnmann transfer w/ simultaneous plane change)
a1 = radius; % m, radius of Ceres
aT = (a1+a2)/2; % m, semimajor axis of transfer ellipse
h = sqrt(2*mu*a1*a2/(a1+a2)); % kg-m2/s, angular momentum of transfer orbit
VA = sqrt(mu/a2); % m/s, velocity of circular orbit at LCO
%VB = sqrt(mu/a1); % m/s, velocity of circular orbit at surface of Ceres
%Vs = 92.61; % m/s, velocity of surface of Ceres at equator
Vap = h/a2; % m/s, velocity at apoapsis
Vpe = h/a1; % m/s, velocity at periapsis
%delVA = VA*(1-sqrt(2*a1/(a1+a2))); % m/s, depart LCO
%delVB = VB*(sqrt(2*a2/(a1+a2))-1); % m/s, capture at surface orbit

% total delta V's
%delVpc = 2*VA*sin(axial/2); % m/s, plane change at LCO to equatorial plane
delVdep = sqrt(Vap^2+VA^2 - 2*Vap*VA*cos(ang)); % m/s, simultaneous plane change and departure burn at LCO
delVcap = Vpe;%-Vs; % m/s, capture
delVtot_dec = delVdep + delVcap;

% time of transfer
T = pi*sqrt(aT^3/mu); % sec, time to transfer
fprintf(' Descent\n')
fprintf('   delta V departure: %.2f m/s\n',delVdep)
fprintf('   delta V capture: %.2f m/s\n',delVcap)
fprintf('   delta V total: %.2f m/s\n',delVtot_dec)
fprintf('   Descent Transfer Time: %.2f hrs\n',T/3600)

% descent mass ratios
Isp = 320; % sec, specific impulse
M0_12 = 2500; % kg, wet mass before descent
Mratio = exp(delVdep/(g*Isp));
Mra = Mratio;
M_D = M0_12/Mratio;
Mratio = exp(delVcap/g/Isp);
Mra = Mra*Mratio;
M_C = M_D/Mratio;

fprintf('   Landing Mass: %.2f kg\n',M_C)
req_prop = M0_12-M_C;
fprintf('   Required Propellant Mass: %.2f kg\n',req_prop)

% surface to low Ceres orbit (hohmann transfer w/ simultaneous plane change)
aT = (a1+a2)/2; % m, semimajor axis of transfer ellipse
h = sqrt(2*mu*a1*a2/(a1+a2)); % kg-m2/s, angular momentum of transfer orbit
%VA = sqrt(mu/a1); % m/s, velocity of circular orbit at surface of Ceres
VB = sqrt(mu/a2); % m/s, velocity of circular orbit at LCO
%Vs = 92.61; % m/s, velocity of surface of Ceres at equator
Vap = h/a2; % m/s, velocity at apoapsis
Vpe = h/a1; % m/s, velocity at periapsis
%delVB = VB*(1-sqrt(2*a1/(a1+a2))); % m/s, capture at LCO

% total delta V's
%delVleave = ; % m/s, surface burn to leave surface and enter surface orbit
%delVpc = 2*VA*sin(axial/2); % m/s, plane change at LCO to elliptical plane
delVdep = Vpe; %- Vs; % m/s, departure
delVcap = sqrt(Vap^2+VB^2-2*Vap*VB*cos(ang)); % m/s, simultaneous plane change and capture burn at LCO
delVtot = delVdep + delVcap;

% time of transfer
T = pi*sqrt(aT^3/mu); % sec, time to transfer

fprintf('\n')
fprintf(' Ascent\n')
fprintf('   delta V departure: %.2f m/s\n',delVdep)
fprintf('   delta V capture: %.2f m/s\n',delVcap)
fprintf('   delta V total: %.2f m/s\n',delVtot)
fprintf('   Descent Transfer Time: %.2f hrs\n',T/3600)

% ascent mass ratios
M_left = 100; % kg, mass left on Ceres
M0_21 = M_C - M_left; % kg, mass before ascent
Mratio = exp(delVdep/g/Isp);
Mra = Mra*Mratio;
M_D = M0_21/Mratio;
Mratio = exp(delVcap/g/Isp);
Mra = Mra*Mratio;
M_C = M_D/Mratio;

fprintf('   Rendezvous Mass: %.2f kg\n',M_C)
fprintf('   Required Propellant Mass: %.2f kg\n',M0_21-M_C)
fprintf('\n')
fprintf(' Overall\n')
fprintf('   Overall delta V total: %.2f m/s\n',delVtot+delVtot_dec)
fprintf('   Total Required Propellant Mass: %.2f kg\n\n',M0_21-M_C + req_prop)

% plot final orbit
craftmass = M0_12;
vel = VB/1000;
orbits_ceres(craftmass,angdeg,vel)


%% General Surface Study Mission
fprintf('\n')
fprintf('General Surface Study Mission: General Lander 1, 2, 3, & 4\n')
fprintf('\n')
% landing site = 30 degrees from equator (both northern and southern hemisphere)
angdeg = 30+axial;
ang = deg2rad(angdeg);

% low Ceres orbit --> surface (hohnmann transfer w/ simultaneous plane change)
a1 = radius; % m, radius of Ceres
aT = (a1+a2)/2; % m, semimajor axis of transfer ellipse
h = sqrt(2*mu*a1*a2/(a1+a2)); % kg-m2/s, angular momentum of transfer orbit
VA = sqrt(mu/a2); % m/s, velocity of circular orbit at LCO
%VB = sqrt(mu/a1); % m/s, velocity of circular orbit at surface of Ceres
Vs = 92.61; % m/s, velocity of surface of Ceres at equator
wc = Vs/radius; % angular velocity of Ceres rotation
Vgr = wc*(radius*cos(ang));
Vap = h/a2; % m/s, velocity at apoapsis
Vpe = h/a1; % m/s, velocity at periapsis
%delVA = VA*(1-sqrt(2*a1/(a1+a2))); % m/s, depart LCO
%delVB = VB*(sqrt(2*a2/(a1+a2))-1); % m/s, capture at surface orbit

% total delta V's
%delVpc = 2*VA*sin(axial/2); % m/s, plane change at LCO to equatorial plane
delVdep = sqrt(Vap^2+VA^2 - 2*Vap*VA*cos(ang)); % m/s, simultaneous plane change and departure burn at LCO
delVcap = Vpe-Vgr; % m/s, capture
delVtot_desc = delVdep + delVcap;

% time of transfer
T = pi*sqrt(aT^3/mu); % sec, time to transfer
fprintf(' Descent\n')
fprintf('   delta V departure: %.2f m/s\n',delVdep)
fprintf('   delta V capture: %.2f m/s\n',delVcap)
fprintf('   delta V total: %.2f m/s\n',delVtot_desc)
fprintf('   Descent Transfer Time: %.2f hrs\n',T/3600)

% descent mass ratios
Isp = 320; % sec, specific impulse
M0_12 = 2500; % kg, wet mass before descent
Mratio = exp(delVdep/(g*Isp));
Mra = Mratio;
M_D = M0_12/Mratio;
Mratio = exp(delVcap/g/Isp);
Mra = Mra*Mratio;
M_C = M_D/Mratio;

fprintf('   Landing Mass: %.2f kg\n',M_C)
req_prop = M0_12-M_C;
fprintf('   Required Propellant Mass: %.2f kg\n',req_prop)

% surface to low Ceres orbit (hohmann transfer w/ simultaneous plane change)
aT = (a1+a2)/2; % m, semimajor axis of transfer ellipse
h = sqrt(2*mu*a1*a2/(a1+a2)); % kg-m2/s, angular momentum of transfer orbit
%VA = sqrt(mu/a1); % m/s, velocity of circular orbit at surface of Ceres
VB = sqrt(mu/a2); % m/s, velocity of circular orbit at LCO
Vs = 92.61; % m/s, velocity of surface of Ceres at equator
wc = Vs/radius; % angular velocity of Ceres rotation
Vgr = wc*(radius*cos(ang));
Vap = h/a2; % m/s, velocity at apoapsis
Vpe = h/a1; % m/s, velocity at periapsis
delVB = VB*(1-sqrt(2*a1/(a1+a2))); % m/s, capture at LCO

% total delta V's
%delVleave = ; % m/s, surface burn to leave surface and enter surface orbit
%delVpc = 2*VA*sin(axial/2); % m/s, plane change at LCO to elliptical plane
delVdep = Vpe-Vgr; % m/s, departure
delVcap = sqrt(Vap^2+VB^2-2*Vap*VB*cos(ang)); % m/s, simultaneous plane change and capture burn at LCO
delVtot = delVdep + delVcap;

% time of transfer
T = pi*sqrt(aT^3/mu); % sec, time to transfer

fprintf('\n')
fprintf(' Ascent\n')
fprintf('   delta V departure: %.2f m/s\n',delVdep)
fprintf('   delta V capture: %.2f m/s\n',delVcap)
fprintf('   delta V total: %.2f m/s\n',delVtot)
fprintf('   Descent Transfer Time: %.2f hrs\n',T/3600)

% ascent mass ratios
M_left = 100; % kg, mass left on Ceres
M0_21 = M_C - M_left; % kg, mass before ascent
Mratio = exp(delVdep/g/Isp);
Mra = Mra*Mratio;
M_D = M0_21/Mratio;
Mratio = exp(delVcap/g/Isp);
Mra = Mra*Mratio;
M_C = M_D/Mratio;

fprintf('   Rendezvous Mass: %.2f kg\n',M_C)
fprintf('   Required Propellant Mass: %.2f kg\n',M0_21-M_C)
fprintf('\n')
fprintf(' Overall\n')
fprintf('   Overall delta V total: %.2f m/s\n',delVtot+delVtot_desc)
fprintf('   Total Required Propellant Mass: %.2f kg\n\n',M0_21-M_C + req_prop)

% plot final orbit
craftmass = M0_12;
vel = VB/1000;
orbits_ceres(craftmass,angdeg,vel) % lander 1
orbits_ceres(craftmass,-angdeg,vel) % lander 2

%% Surface Mapping Study Mission
fprintf('\n')
fprintf('Surface Mapping Study Mission: Surface Mapping Orbiter\n')
fprintf('\n')
% desired orbit = polar orbit
angdeg = 90+axial;
ang = deg2rad(angdeg);

VA = sqrt(mu/a2); % m/s, velocity of circular orbit at LCO

% total delta V's
delVpc = 2*VA*sin(ang/2); % m/s, plane change at LCO
delVdep = delVpc;

% time of transfer
T = pi*sqrt(aT^3/mu); % sec, time to transfer
fprintf(' Plane Change\n')
fprintf('   delta V total: %.2f m/s\n',delVdep)
fprintf('   Plane Change Time: %.2f hrs\n',T/3600)

% descent mass ratios
Isp = 3000; % sec, specific impulse
M0_12 = 1500; % kg, wet mass before descent
Mratio = exp(delVdep/(g*Isp));
Mra = Mratio;
M_D = M0_12/Mratio;
Mratio = exp(delVcap/g/Isp);
Mra = Mra*Mratio;
M_C = M_D/Mratio;

fprintf('   Final Mass: %.2f kg\n',M_C)
fprintf('   Required Propellant Mass: %.2f kg\n\n',M0_12-M_C)

% plot final orbit
craftmass = M0_12;
vel = VB/1000;
orbits_ceres(craftmass,angdeg,vel)

%% Communications Array Mission
fprintf('\n')
fprintf('Communications Array Mission: Comms Orbiter 1 & Comms Orbiter 2\n')
fprintf('\n')

% desired orbit: high elliptic orbit (HEO) with 60 deg axial inclination
angdeg = 30+axial;
ang = deg2rad(angdeg);

% low Ceres orbit --> HEO with 60 deg axial inclination
a1 = radius; % m, radius of Ceres
rA = soi*(0.03); % radius of apoapsis
rP = a2; % radius of periapsis
aT = (rA+rP)/2; % m, semimajor axis of transfer ellipse
e = (rA-rP)/(rA+rP); % eccentricity
h = sqrt(2*mu*rA*rP/(rA+rP)); % kg-m2/s, angular momentum of HEO
VA = sqrt(mu/rP); % m/s, velocity of circular orbit at LCO
Vap = h/rA; % m/s, velocity at apoapsis
Vpe = h/rP; % m/s, velocity at periapsis

% total delta V's
delVdep = sqrt(Vap^2+VA^2 - 2*Vap*VA*cos(ang)); % m/s, simultaneous plane change and orbit change

% time of transfer
T = pi*sqrt(aT^3/mu); % sec, time to transfer
fprintf(' Orbit and Plane Change\n')
fprintf('   delta V total: %.2f m/s\n',delVdep)
fprintf('   Transfer Time: %.2f hrs\n',T/3600)

% descent mass ratios
Isp = 3000; % sec, specific impulse
M0_12 = 1500; % kg, wet mass before orbit and place changes
Mratio = exp(delVdep/(g*Isp));
Mra = Mratio;
M_D = M0_12/Mratio;
Mratio = exp(delVcap/g/Isp);
Mra = Mra*Mratio;
M_C = M_D/Mratio;

fprintf('   Landing Mass: %.2f kg\n',M_C)
fprintf('   Required Propellant Mass: %.2f kg\n',M0_12-M_C)

% plot final orbit
craftmass = M0_12;
vel = Vpe/1000;
orbits_ceres(craftmass,angdeg,vel) % comms orbiter 1
orbits_ceres(craftmass,-angdeg,vel) % comms orbiter 2
orbits_ceres(craftmass,180+angdeg,vel) % comms orbiter 3
orbits_ceres(craftmass,180-angdeg,vel) % comms orbiter 4