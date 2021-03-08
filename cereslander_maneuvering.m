%% Ceres Lander Maneuvering Calculations
clear; clc; close all;

% parameters
mu = 62630000000; % m3/s2, ceres gravitation parameter
axial = 0.06981317008; % rad, axial tilt (or what ever angle from the elliptical plane you want)
radius = 469730; % m, average radius of surface
g = 9.807; % m/s2, acceleration of Earth

% low ceres orbit (change altitude to determine best height)
alt = 100000; % altitude of LCO
a2 = radius + alt; % m, radius of circular orbit at surface

%% low ceres orbit --> surface (hohnmann transfer)
a1 = radius; % m, radius of circular orbit at surface
aT = (a1+a2)/2; % m, semimajor axis of transfer ellipse
h = sqrt(2*mu*a1*a2/(a1+a2)); % kg-m2/s, angular momentum of transfer orbit
VA = sqrt(mu/a2); % m/s, velocity of circular orbit at LCO
VB = sqrt(mu/a1); % m/s, velocity of circular orbit at surface of ceres
Vs = 92.61; % m/s, velocity of surface of ceres at equator
Vap = h/a2; % m/s, velocity at apoapsis
Vpe = h/a1; % m/s, velocity at periapsis
delVA = VA*(1-sqrt(2*a1/(a1+a2))); % m/s, depart LCO
delVB = VB*(sqrt(2*a2/(a1+a2))-1); % m/s, capture at surface orbit


% total delta V's
delVmatch = VB-Vs; % m/s, change velocity to match surface of ceres
delVpc = 2*VA*sin(axial/2); % m/s, plane change at LCO to equatorial plane
delVdep = sqrt(Vap^2+VA^2 - 2*Vap*VA*cos(axial)); % m/s, simultaneous plane change and departure burn at LCO
delVcap = delVB + delVmatch; % m/s, capture at surface orbit
delVtot = delVdep + delVcap;

% time of transfer
T = pi*sqrt(aT^3/mu); % sec, time to transfer

display(delVtot)
display(T)

% descent mass ratios
Isp = 320; % sec, specific impulse
M0_12 = 1966.81; % kg, wet mass before descent
Mratio = exp(delVdep/(g*Isp));
Mra = Mratio;
M_D = M0_12/Mratio;
Mratio = exp(delVcap/g/Isp);
Mra = Mra*Mratio;
M_C = M_D/Mratio;


%% surface to low ceres orbit (hohmann transfer)
a2 = 569730; % m, radius of circular orbit at LCO
aT = (a1+a2)/2; % m, semimajor axis of transfer ellipse
h = sqrt(2*mu*a1*a2/(a1+a2)); % kg-m2/s, angular momentum of transfer orbit
VA = sqrt(mu/a1); % m/s, velocity of circular orbit at surface of ceres
VB = sqrt(mu/a2); % m/s, velocity of circular orbit at LCO
Vs = 92.61; % m/s, velocity of surface of ceres at equator
Vap = h/a2; % m/s, velocity at apoapsis
Vpe = h/a1; % m/s, velocity at periapsis
delVA = VA*(sqrt(2*a2/(a1+a2))-1); % m/s, depart surface orbit
delVB = VB*(1-sqrt(2*a1/(a1+a2))); % m/s, capture at LCO

% total delta V's
delVleave = VA - Vs; % m/s, surface burn to leave surface and enter surface orbit
delVpc = 2*VA*sin(axial/2); % m/s, plane change at LCO to elliptical plane
delVdep = delVleave + delVA; % m/s, enter transfer orbit from surface orbit
delVcap = sqrt(Vap^2+VB^2-2*Vap*VB*cos(axial)); % m/s, simultaneous plane change and capture burn at LCO
delVtot = delVdep + delVcap;

% time of transfer
T = pi*sqrt(aT^3/mu); % sec, time to transfer

display(delVtot)
display(T)

% ascent mass ratios
M_left = 100; % kg, mass left on ceres
M0_21 = M_C - M_left; % kg, wet mass before ascent
Mratio = exp(delVdep/g/Isp);
Mra = Mra*Mratio;
M_D = M0_21/Mratio;
Mratio = exp(delVcap/g/Isp);
Mra = Mra*Mratio;
M_C = M_D/Mratio;

display(M_C)
display(Mra)
