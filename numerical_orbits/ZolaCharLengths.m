%% Zola Characteristic Lengths
% 5.2.21
% sdf33
clear
clc
close all

%...Project Demeter
%...low-thrust trajectory from Earth to Ceres, and return
%...constant thrust
%...note: t_ , time, values start from t = 0

%...General Parameters
g = 9.807; % m/s2, Earth gravity

%...Alpha Parameters
mo = 325000; % kg, initial wet mass
mp = 300000; % kg, propellant mass
F = 102.1; % N, thrust
Isp = 8180; % sec, specific impulse
vj = Isp*g; % effective velocity


%%...Phase 1...accelerating
t1 = 100; % days, time
t1 = t1*24*3600; % secs, time
tphase1 = 0:1:t1;
a1_hist = ao / (1 - ao*tphase1/vj); % m/s2, acceleration history as a function of t, for 0 < t < t1
a1 = ao / (1 - ao*t1/vj); % m/s2, acceleration magnitude, for 0 < t < t1
V1_hist = -vj*ln(1-ao*tphase1/vj); % m/s, velocity history as a function of t, for 0 < t < t1
V1 = -vj*ln(1-ao*t1/vj); % m/s, velocity magnitude, for 0 < t < t1
L1 = (vj^2/ao)*((1-ao*t1/vj)*ln(1-ao*t1/vj) + ao*t1/vj); % m, characteristic length of Phase 1

%%...Phase 2...coasting
t2 = 400; % days, time
t2 = t2*24*3600; % secs, time
tphase2 = t1:1:t1+t2;
Vmax = V1; % m/s, maximum velocity during coasting, for t1 < t < t2
L2 = Vmax*(t2-t1); % m, characteristic length, for t1 < t < t2

%%...Phase 3...deccelerating
T = 100; % days, time
T = T*24*3600; % secs, time
tphase3 = t2:1:T;
a3_hist = ao/(1 - ao*t1/vj - ao*(t-t2)/vj); % m/s2, accleration history as a function of t, for t2 < t < T
V3_hist = -vj*ln(1-ao*tphase3/vj); % m/s, velocity history as a function of t, for 0 < t < t1
L3 = vj^2/ao*((ai*t1/vj)^2 - (1-ao*t1/vj)*ln(1-ao*t1/vj) - ao*t1/vj); % m, characteristic length, for t2 < t < T

%%...Mission Totals
L = L1+L2+L3;


figure
hold on
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
plot(tphase1,V1_hist)
plot(tphase2,Vmax)
plot(tphase3,)





