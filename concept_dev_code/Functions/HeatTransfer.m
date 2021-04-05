function [Q_t] = HeatTransfer(D,T,R,A)
%% Note
% D is the is the distance between spacecraft and the planet/star [AU]
% T is the effective temperature of theplanet or the surface temperature of
% the star.
% R is the radius of the planet/star [km]
% A is the Surface Area of the Object
%% Constant
sigma = 5.7e-8; % Stefan-Bolzmann Constant (W/m^2/K^4)
a_u = 149597870691; % AU to m
%% Unit Conversion
D = D.*a_u; % AU -> m
R = R.*1000; % km -> m
%% Heat Flux
G = sigma*T^4.*(R./D);
Q_t = G.*A;
end