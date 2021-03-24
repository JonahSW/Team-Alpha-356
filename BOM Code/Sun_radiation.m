function [distance_from_sun, solar_flux] = Sun_radiation(a_1,a_2)
%Input
% a_1 -> End distance away from the sun (AU)
% a_2 -> End distance away from the sun (AU)
%Output
% Solar_flux of the sun -> W/m^2


%Solar Constants
T_sun = 5772; %temperature of the sun in K
R_sun = 6.957e8; %Solar Equatorial Radius (m)
sigma = 5.67037e-8;% Stefan-Boltzmann Constant (W/m^2*K^4)
solar_surface_flux = sigma*T_sun^4; % Solar Radiation Flux at surface W/m2
solar_constant = 1.3608e3; % Solar Constant (W/m2) <- solar radiation flux at 1 AU

AU = 1.496e11;%AU in m
if a_1 == a_2
    distance_from_sun = a_1*AU*ones(2650);
else
    distance_from_sun = a_1*AU:(((a_2*AU)-(a_1*AU))/2650):a_2*AU;
end 
solar_flux = [];
for i = 1:length(distance_from_sun)
    solar_flux(i) = solar_surface_flux*(R_sun/distance_from_sun(i))^2;
end 
end 