%mdp81, jjs280
%03/21/2021
%Code for calculating solar radiation flux as function of distance from the sun

%% Calculating Temp as a function of Distance from the Sun 
function solar_flux = solar_heat_flux(distance_from_sun)
%Solar Constants
T_sun = 5772; %temperature of the sun in K
R_sun = 6.957e8; %Solar Equatorial Radius (m)
sigma = 5.67037e-8;% Stefan-Boltzmann Constant (W/m^2*K^4)
solar_surface_flux = sigma*T_sun^4; % Solar Radiation Flux at surface W/m2
solar_constant = 1.3608e3; % Solar Constant (W/m2) <- solar radiation flux at 1 AU

AU = 1.496e11;%AU in m

solar_flux = solar_surface_flux*(R_sun./distance_from_sun).^2;

%% Plot Results
plotme = 1;
if plotme == 1
    figure(1)
    semilogy(distance_from_sun,solar_flux);
    xlabel('Distance From Sun (m)')
    ylabel('Solar Flux (w/m^2)')
    title('Solar Radiation Flux vs. Distance From Sun')
    xline(AU,'-','1 AU');
    yline(solar_constant,'-','Solar Constant');
end

end