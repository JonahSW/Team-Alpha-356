%% Zola Characteristic Lengths
% 5.2.21
% sdf33
clear
clc
close all
%...Project Demeter
%...constant thrust

% General Parameters
g = 9.807; % m/s2, Earth gravity

%...Alpha Parameters
mo = 325000; % kg, initial wet mass
mp = 300000; % kg, propellant mass
T = 102.1; % N, thrust
Isp = 8180; % sec, specific impulse
vj = Isp*g; % effective velocity

a = ao / (1 - ao*ta/vj); % m/s2, acceleration as a function of time