%jjs280
%03/17/2021
%Implements a function that can calculate the mass ratio and delta V for a combined orbital plane
%change and elliptic transfer

function [Mratio, deltaV] = elliptic_transfer_with_plane_change(theta,V,Isp)
    
    %Internally Defined Constants
    g0 = 9.807;%Earth Average gravity [m/s^2]
    
    deltaV = 2*V*sin(theta/2);
    Mratio = e^(deltaV/(g0/Isp));
end