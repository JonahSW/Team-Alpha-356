%jjs280
%03/12/2021
%Implements a function that can calculate the travel time, mass ratios, and delta Vs for an elliptic transfer orbit

%a1, a2, are semi major axes of initial and final orbits
%thetaB_prime is the true anomaly of the target planet on arrival
%Isp is the specific impulse of the spacecraft
function [Mratio_A, Mratio_B, deltaV_A, deltaV_B, tT] = elliptic_transfer(a1, a2, thetaB_prime, Isp)
    
    %Internally Defined Constants
    g0 = 9.807;%Earth Average gravity [m/s^2]
    mu_sun = 1.32712e11;%Solar Gravitational Parameter [km^3/sec^2]
    
    V_A = sqrt(mu_sun/a1);%Tangential velocity of initial orbit (circular)
    V_B = sqrt(mu_sun/a2);%Tangential velocity of final orbit (circular)
    
    %Calculations
    eT = (a2-a1)/(a1-a1*cos(thetaB_prime));%Eccentricity of transfer orbit
    V_B_prime = V_B*sqrt(2-(a1/a2)*(1-eT));%Velocity when transfer orbit crosses target orbit
    phi = asin((a1/a2)*sqrt(((1+eT)/(1-eT))/((a1/a2)*(2/(1-eT))-1)));%Flight path angle
    E = 2*atan(sqrt((1-eT)/(1+eT))*tan(thetaB_prime/2));%Eccentric anomaly
    M = E-eT*sin(E);%Mean Anomaly
    
    %Results
    deltaV_A = V_A*(sqrt(1+eT)-1);%departure burn deltaV 
    deltaV_B = sqrt(V_B^2 + V_B_prime^2 - 2*V_B*V_B_prime*cos(phi - (pi/2)));
    Mratio_A = exp(deltaV_A/(g0*Isp));
    Mratio_B = exp(deltaV_B/(g0*Isp));
    tT = M*sqrt(a1^3/(mu_sun*(1-eT)^3))*(1/(24*3600));
    
end