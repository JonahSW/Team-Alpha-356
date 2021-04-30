%jjs280
%03/31/2021
%Wrapper function for testing cowell's method code

positions = planetary_positions(1);


%Implements a function that calculates the wait time between hohmann transfers

%a1, a2, are semi major axes of initial and final orbits
function [t_wait] = wait_time(a1, a2)
    mu_sun = 1.32712e20;% Solar Gravitational Parameter [m^3/sec^2]
    %Convert to m from AU
    AU = 1.496e11;%AU in m
    a1 = a1*AU;
    a2 = a2*AU;
    %Calculations
    alpha = pi*(((a1+a2)/(2*a1))^(3/2)-1);
    S = (((2*pi)*sqrt(a1^3/mu_sun))^(-1) - ((2*pi)*sqrt(a2^3/mu_sun))^(-1))^(-1);%Synodic Period
    t_wait = (2-(alpha/pi))*(S/(24*3600));
end