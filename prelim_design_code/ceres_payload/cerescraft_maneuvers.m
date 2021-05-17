
%% Ceres Craft Maneuvers
clear; clc; close all;

% sdf33
% updated 5.5.21

% select initial altitude
% assume all spacecrafts begin at LCO in polar orbit
% finds delta V, mass, mass ratios, duration for each transfer

%% Missions
% Polar Exploration Mission
% Surface Exploration Mission
% Surface Mapping Mission
% Communications Array Mission

%% Initial Parameters
% parameters
mu = 62630000000; % m3/s2, Ceres gravitation parameter
axial = 0.06981317008; % rad, axial tilt
radius = 469730; % m, average radius of surface
g = 9.807; % m/s2, acceleration of Earth
soi = 7.7e7; % m

% low Ceres orbit
alt = 700000; % m, altitude of LCO
a2 = radius + alt; % m, radius of LCO

%% Polar Exploration Mission
fprintf('Polar Exploration Mission\n\n')

% landing site = north pole
angdeg = 20+axial;
ang = deg2rad(angdeg);

% low Ceres orbit --> surface (hohnmann transfer w/ simultaneous plane change)
a1 = radius; % m, radius of Ceres
aT = (a1+a2)/2; % m, semimajor axis of transfer ellipse
h = sqrt(2*mu*a1*a2/(a1+a2)); % kg-m2/s, angular momentum of transfer orbit
VA = sqrt(mu/a2); % m/s, velocity of circular orbit at LCO
Vap = h/a2; % m/s, velocity at apoapsis
Vpe = h/a1; % m/s, velocity at periapsis

% total delta V's
delVdep = sqrt(Vap^2+VA^2 - 2*Vap*VA*cos(ang)); % m/s, simultaneous plane change and departure burn at LCO
delVcap = Vpe;%-Vs; % m/s, capture
delVtot_dec = delVdep + delVcap;

% time of orbit change
T = pi*sqrt(aT^3/mu); % sec, time to transfer
fprintf(' Descent\n')
fprintf('   delta V departure: %.2f m/s\n',delVdep)
fprintf('   delta V capture: %.2f m/s\n',delVcap)
fprintf('   delta V total: %.2f m/s\n',delVtot_dec)
fprintf('   Descent Transfer Time: %.2f hrs\n',T/3600)

% descent mass ratios
Isp = 321; % sec, specific impulse
M0_12 = 1500; % kg, wet mass before descent
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
VB = sqrt(mu/a2); % m/s, velocity of circular orbit at LCO
Vap = h/a2; % m/s, velocity at apoapsis
Vpe = h/a1; % m/s, velocity at periapsis

% total delta V's
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
M_left = 161; % kg, mass left on Ceres
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


%% Equatorial Exploration Mission
fprintf('\n')
fprintf('Equatorial Exploration Mission\n')
fprintf('\n')
% landing site = 30 degrees from equator (both northern and southern hemisphere)
angdeg = 80+axial;
ang = deg2rad(angdeg);

% low Ceres orbit --> surface (hohnmann transfer w/ simultaneous plane change)
a1 = radius; % m, radius of Ceres
aT = (a1+a2)/2; % m, semimajor axis of transfer ellipse
h = sqrt(2*mu*a1*a2/(a1+a2)); % kg-m2/s, angular momentum of transfer orbit
VA = sqrt(mu/a2); % m/s, velocity of circular orbit at LCO
Vs = 92.61; % m/s, velocity of surface of Ceres at equator
wc = Vs/radius; % angular velocity of Ceres rotation
Vgr = wc*(radius*cos(ang));
Vap = h/a2; % m/s, velocity at apoapsis
Vpe = h/a1; % m/s, velocity at periapsis

% total delta V's
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
M0_12 = 1800; % kg, wet mass before descent
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
VB = sqrt(mu/a2); % m/s, velocity of circular orbit at LCO
Vs = 92.61; % m/s, velocity of surface of Ceres at equator
wc = Vs/radius; % angular velocity of Ceres rotation
Vgr = wc*(radius*cos(ang));
Vap = h/a2; % m/s, velocity at apoapsis
Vpe = h/a1; % m/s, velocity at periapsis
delVB = VB*(1-sqrt(2*a1/(a1+a2))); % m/s, capture at LCO

% total delta V's
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
M_left = 161; % kg, mass left on Ceres
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



%% Surface Mapping Mission
fprintf('\n')
fprintf('Surface Mapping Mission\n')
fprintf('\n')

% low Ceres orbit --> surface mapping orbit (hohnmann transfer)
a1 = 100000+radius; % m, surface mapping orbit radius
aT = (a1+a2)/2; % m, semimajor axis of transfer ellipse
h = sqrt(2*mu*a1*a2/(a1+a2)); % kg-m2/s, angular momentum of transfer orbit
VA = sqrt(mu/a2); % m/s, velocity of circular orbit at LCO
VB = sqrt(mu/a1); % m/s, velocity of circular orbit at surface mapping orbit
Vap = h/a2; % m/s, velocity at apoapsis
Vpe = h/a1; % m/s, velocity at periapsis

% total delta V's
delVdep = VA*(1 - sqrt(2*a1/(a1+a2))); % m/s, departure
delVcap = VB*(sqrt(2*a2/(a1+a2)) - 1); % m/s, capture
delVtot = delVdep + delVcap; % m/s, total

% time of transfer
T = pi*sqrt(aT^3/mu); % sec, time to transfer
fprintf(' Plane Change\n')
fprintf('   delta V total: %.2f m/s\n',delVtot)
fprintf('   Plane Change Time: %.2f hrs\n',T/3600)

% descent mass ratios
Isp = 321; % sec, specific impulse
M0_12 = 1000; % kg, wet mass before descent
Mratio = exp(delVdep/(g*Isp));
Mra = Mratio;
M_D = M0_12/Mratio;
Mratio = exp(delVcap/g/Isp);
Mra = Mra*Mratio;
M_C = M_D/Mratio;

fprintf('   Final Mass: %.2f kg\n',M_C)
fprintf('   Required Propellant Mass: %.2f kg\n\n',M0_12-M_C)


%% Communications Mission
fprintf('\n')
fprintf('Communications Mission\n')
fprintf('\n')

% desired orbit: highly elliptical orbit (HEO)
angdeg = 30+axial;
ang = deg2rad(angdeg);

% low Ceres orbit --> HEO with 60 deg axial inclination
a1 = radius; % m, radius of Ceres
rA = 2000000; % m, radius of apoapsis
rP = a2; % m, radius of periapsis
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
Isp = 321; % sec, specific impulse
M0_12 = 30; % kg, wet mass before orbit and place changes
Mratio = exp(delVdep/(g*Isp));
Mra = Mratio;
M_D = M0_12/Mratio;
Mratio = exp(delVcap/g/Isp);
Mra = Mra*Mratio;
M_C = M_D/Mratio;

fprintf('   Final Mass: %.2f kg\n',M_C)
fprintf('   Required Propellant Mass: %.2f kg\n',M0_12-M_C)

