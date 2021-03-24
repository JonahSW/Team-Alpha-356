function [T] = EffectiveTemp(a_1,Ba)
%% Note
% D is the distance from the Planet to the Sun (AU);
% Ba is the Bond albedo of the Planet
% Ref. http://en.wikipedia.org/wiki/Bond_albedo
%% Constant
sigma = 5.7e-8; % Stefan-Bolzmann Constant (W/m^2/K^4)
L_sun = 3.827E26; % Sun Luminosity (W)
a_u = 149597870691; % AU to m
%% Effective Temperature of the Planet due to Sun
T = ( L_sun*(1-Ba) / (16*pi*sigma*(a_1*a_u)^2) )^(1/4);
end