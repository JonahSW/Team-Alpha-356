%% Zola Characteristic Lengths
% 5.2.21
% sdf33
clear; clc; close all

%...Project Demeter
%...low-thrust trajectory from Earth to Ceres, and return
%...constant thrust
%...note: t_ , time, values start from t = 0

%...general parameters
g = 9.807; % m/s2, Earth gravity

%...Alpha parameters
mo = 325000; % kg, initial wet mass
mp = 300000; % kg, propellant mass
F = 102.1; % N, thrust
ao = F/mo; % m/s2, initial acceleration
Isp = 8180; % sec, specific impulse
vj = Isp*g; % effective velocity


%%...phase 1...accelerating, 0 < t < t1
t1 = 100; % days, time
t1 = t1*24*3600; % secs, time
tphase1 = 0:1:t1;
a1_hist = ao./(1 - ao.*tphase1./vj); % m/s2, acceleration history as a function of t
a1 = ao./(1 - ao.*t1./vj); % m/s2, acceleration magnitude
V1_hist = -vj.*log(1-ao.*tphase1./vj); % m/s, velocity history as a function of t
V1 = -vj.*log(1-ao.*t1./vj); % m/s, velocity magnitude
L1 = (vj.^2./ao).*((1-ao.*t1./vj).*log(1-ao.*t1./vj) + ao.*t1./vj); % m, characteristic length

%%...phase 2...coasting, t1 < t < t2
t2 = t1+200; % days, time
t2 = t2*24*3600; % secs, time
tphase2 = t1:1:t2;
Vmax = V1; % m/s, maximum velocity during coasting
L2 = Vmax.*(t2-t1); % m, characteristic length

%%...phase 3...deccelerating, t2 < t < T
T = t2+100; % days, time
T = T*24*3600; % secs, time
tphase3 = t2:1:T;
a3_hist = ao./(1-ao.*t1./vj-ao.*(tphase3-t2)./vj); % m/s2, accleration history as a function of t
V3_hist = -vj.*log(1-ao.*tphase3./vj); % m/s, velocity history as a function of t
L3 = vj.^2./ao.*((ao.*t1./vj).^2-(1-ao.*t1./vj).*log(1-ao.*t1./vj) - ao.*t1./vj); % m, characteristic length

%%...mission totals
L = L1+L2+L3;

%%...plotting
figure
hold on
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
plot(tphase1,V1_hist)
plot(tphase2,Vmax)
plot(tphase3,V3_hist)





