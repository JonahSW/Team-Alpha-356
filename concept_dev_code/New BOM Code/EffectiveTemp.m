function [T_2] = EffectiveTemp(a_1,a_2,a_earth,a_ceres,distance_from_sun,surface_area,T_1)
%% Note
% D is the distance from the Planet to the Sun (AU);
% Ba is the Bond albedo of the Planet
% Ref. http://en.wikipedia.org/wiki/Bond_albedo
%% Constant

%Ba = 0.0331;
Ba_earth = 0.3;
Ba_ceres = 0.09;
E = 0.0331;
R_sun = 6.957e8;

sigma = 5.7e-8; % Stefan-Bolzmann Constant (W/m^2/K^4)
T_sun = 5772; %temperature of the sun in K
solar_surface_flux = sigma*T_sun^4;
solar_flux = solar_surface_flux.*(R_sun./distance_from_sun).^2;
Temp_2 = ((solar_flux./(E*sigma*surface_area)) + (T_1^4)).^(1/4);
L_sun = 3.827E26; % Sun Luminosity (W)

%% Effective Temperature of the Planet due to Sun

if a_1 == a_2 && a_1 == a_earth
    T_2 = Temp_2 + (( L_sun*(1-Ba_earth)./ (16*pi*sigma*(distance_from_sun).^2)).^(1/4));
elseif a_1 == a_2 && a_1 == a_ceres
    T_2 = Temp_2 + (( L_sun*(1-Ba_ceres)./ (16*pi*sigma*(distance_from_sun).^2)).^(1/4));
else
    T_2 = Temp_2;
end 
end 