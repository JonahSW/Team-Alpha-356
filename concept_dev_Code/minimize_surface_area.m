%jjs280
%03/22/2021
%Function for optimizing the surface area of a cylindrical tank with hemispherical endcaps.
%Incomplete and currently not in use

function [r,h,minimum_SA] = minimize_surface_area(volume)
    %V = (4/3)*pi*r^3 + pi*r^2*h
    %SA = 2*pi*r*h + 4*pi*r^2 = 2V/r - (4/3)*pi*r^2 + 4*pi*r^2
    %dSA/dr = -2*V*r^-2 - (16/3)*pi*r
    %dSA/dr = 0 -> r = (2*V/(-16*pi/3))^(1/3)
    r = (2.*volume/(8.*pi./3)).^(1/3);
    h = (volume - (4./3).*pi.*r.^3)/(pi.*r.^2);
    minimum_SA = 2.*volume./r - (4/3).*pi.*r.^2 + 4.*pi.*r.^2;   
end