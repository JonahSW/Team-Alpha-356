function [t_trans,delV_tot,delV_A,delV_B] = Hohmann(a_1,a_2)
%% Note
% a_1 is the distance from departure planet to sun (AU);
% a_2 is the distance from target planet to sun (AU);
%% Constant
mu_sun = 1.33E+11; % Gravitation Parameter Sun (km^3/s^2)
a_u = 149597870.691; % AU to Km

%% Hohmann Transfer
V_A = sqrt(mu_sun/(a_1 * a_u)); % Va (km/s)
V_B = sqrt(mu_sun/(a_2 * a_u)); % Vb (km/s)
delV_A = abs(V_A * (((2*a_2/(a_1+a_2))^(1/2))-1)); % First Burn (km/s)
delV_B = abs(V_B * (1-((2*a_1/(a_1+a_2))^(1/2)))); % First Burn (km/s)
delV_tot = abs(delV_A + delV_B); % Total Burn Velocity (km/s)

a_T = (a_1+a_2)/2; % Semi Major Axis
t_trans = 0.5*sqrt(a_T^3); % Transfer Time (years)
end