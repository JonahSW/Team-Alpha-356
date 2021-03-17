%jjs280
%03/17/2021
%Implements a function that calculates the duration of a hohmann trajectory

%a1, a2, are semi major axes of initial and final orbits
function [duration] = hohmann_duration(a1, a2)
    mu_sun = 1.32712e20;% Solar Gravitational Parameter [m^3/sec^2]
    aT = (a1+a2)/2;
    duration = pi*sqrt(aT/mu_sun);
    
end