%jjs280
%03/17/2021
%Implements a function that can calculate the mass ratio, delta V, and departure angle psi for 
%a Hohmann Transfer departure burn assuming circular coplanar orbits and using patched conics

%a1, a2, are semi major axes of initial and final orbits
%delta is the plane change desired for the burn
%altitude is the altitude of the starting orbit
%Isp is the specific impulse of the spacecraft
function [Mratio, deltaV_departure, psi] = hohmann_departure_burn(a1, a2, mu_1, delta, altitude, Isp)
    
    %Internally Defined Constants
    g0 = 9.807;% Earth Average surface gravity [m/s^2]
    mu_sun = 1.32712e20;% Solar Gravitational Parameter [m^3/sec^2]

    %Convert to m from AU
    AU = 1.496e11;%AU in m
    a1 = a1*AU;
    a2 = a2*AU;
    
    %Calculations
    %Hohmann deltaV
    V_1 = sqrt(mu_sun/a1);%Tangential velocity of planet A orbit (circular)
    deltaV_1 = V_1*(sqrt((2*a2)/(a2+a1))-1);
    %Patched Conics
    V_start = sqrt(mu_1/altitude);%Tangential velocity of spacecraft orbit (circular)
    a_hyperbola = mu_1/(deltaV_1^2);
    e_hyperbola = (altitude/a_hyperbola)+1;
    b_hyperbola = a_hyperbola*sqrt(e_hyperbola^2-1);
    V_hyperbola = deltaV_1*sqrt((e_hyperbola+1)/(e_hyperbola-1));%Hyperbolic excess velocity
    
    %Accounting for plane change of delta
    deltaV_departure = sqrt(V_hyperbola^2 + V_start^2 - 2*V_hyperbola*V_start*cos(delta));
    Mratio = exp(deltaV_departure/(g0*Isp));
    psi = atan(b_hyperbola/a_hyperbola);

end