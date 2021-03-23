function [A] = TankDimensions(r,n,M_p)
%% Notes
% r is the height of the tank (m);
% n is the number of the tank;
% M_p is the mass of the propplant (MT)

%% Property of Propplant
% Try to use refprop!!!!!
density = 71; % Density of Liquid Hydrogen (kg/m^3)
v = (137E6)/10; % viscosity of Liquid Hydrogen gram/cm sec 10^6 -> kg/m sec
C_p = 9.58*1000; % J/kg*k
kv = v/density; % Kinematic Viscosity of Liquid Hydrogen
%% Unit Conversion
M_p = M_p * 1000; % MT -> kg

%% Tank Dimensions
volume = M_p./density./n;
h = volume./(pi.*r.^2);
A = n*(2*pi.*r.*h)+(2.*pi.*r.^2); % Outter Surface Area of All Tanks

end