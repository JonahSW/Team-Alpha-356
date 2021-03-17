%jjs280
%03/17/2021
%Implements a function that can calculate the mass ratio and delta V for a Hohmann Transfer
%capture burn assuming circular coplanar orbits and using patched conics

%a1, a2, are semi major axes of initial and final orbits
%theta is the plane change desired for the burn
%altitude is the altitude of the final orbit
%Isp is the specific impulse of the spacecraft
function [Mratio, deltaV_capture, psi] = hohmann_capture_burn(a1, a2, mu_2, delta, altitude, Isp)
    
    %Internally Defined Constants
    g0 = 9.807;% Earth Average surface gravity [m/s^2]
    mu_sun = 1.32712e20;% Solar Gravitational Parameter [m^3/sec^2]

    %Convert to m from AU
    AU = 1.496e11;%AU in m
    a1 = a1*AU;
    a2 = a2*AU;
    
    %Calculations
    %Hohmann deltaV
    V_2 = sqrt(mu_sun/a2);%Tangential velocity of planet 2 orbit (circular)
    deltaV_2 = V_2*(1-sqrt((2*a1)/(a2+a1)));
    %Patched Conics
    V_end = sqrt(mu_2/altitude);%Tangential velocity of spacecraft orbit (circular)
    a_hyperbola = mu_2/(deltaV_2^2);
    e_hyperbola = (altitude/a_hyperbola)+1;
    b_hyperbola = a_hyperbola*sqrt(e_hyperbola^2-1);
    V_hyperbola = deltaV_2*sqrt((e_hyperbola+1)/(e_hyperbola-1));%Hyperbolic excess velocity
    
    %Accounting for plane change of delta
    deltaV_capture = sqrt(V_hyperbola^2 + V_end^2 - 2*V_hyperbola*V_end*cos(delta));
    Mratio = exp(deltaV_capture/(g0*Isp));
    psi = atan(b_hyperbola/a_hyperbola);

end