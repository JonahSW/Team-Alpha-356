%jjs280
%03/12/2021
%Implements a function that can calculate the mass ratio and delta V for an orbital plane change

function [Mratio, deltaV] = plane_change(theta,V,Isp)
    
    %Internally Defined Constants
    g0 = 9.807;%Earth Average gravity [m/s^2]
    
    deltaV = 2*V*sin(theta/2);
    Mratio = exp(deltaV/(g0*Isp));
end