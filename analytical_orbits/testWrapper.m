%jjs280
%03/14/2021
%Wrapper file for testing analytical orbit calculation tools
close all
clc
clear

%Load ephemeris data
ephemerides = ephemerides;

a_Earth = 1.496e+8*mean(ephemerides(4,:));%Convert AU to km
a_Ceres = 1.496e+8*mean(ephemerides(8,:));%Convert AU to km
    
Isp = 950;

mu_sun = 1.32712e11;%Solar Gravitational Parameter [km^3/sec^2]
L12_Hohmann = pi*(1-((a_Earth+a_Ceres)/(2*a_Ceres))^(3/2));%Determines starting separation angle for Hohmann transfer
thetaB_prime = pi;%Limiting case Hohmann Transfer

[Mratio_A, Mratio_B, deltaV_A, deltaV_B, tT] = elliptic_transfer(a_Earth, a_Ceres, thetaB_prime, Isp);