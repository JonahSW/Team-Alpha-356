%jjs280
%03/17/2021
%Implements a function that calculates the wait time between hohmann transfers

%a1, a2, are semi major axes of initial and final orbits
function [t_wait] = wait_time(a1, a2, mu1, mu2)
    mu_sun = 1.32712e20;% Solar Gravitational Parameter [m^3/sec^2]
    alpha = pi*(((a1+a2)/(2*a1))^(3/2)-1);
    S = (sqrt(mu1/a1)/(2*pi) - sqrt(mu2/a2)/(2*pi))^-1;%Synodic Period
    t_wait = (1-(alpha/pi))*S;
end