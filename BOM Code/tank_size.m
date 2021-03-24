function [M_w_first, volume_tank_first, h_tank_first, surface_area, thickness, volume_tank, tank_weight] = tank_size(total_mass_fraction, M_d, r_tank, number_of_tanks)
%Input

%Output

%% Calculating Tank Size (Cylinder)
M_w_first = total_mass_fraction*M_d; %Wet Mass [MT]
M_p_first = M_w_first - M_d;
M_P_first_kg = M_p_first *1000; 
liquid_hydrogen = 71; % kg/m^3
volume_tank_first = (M_P_first_kg/liquid_hydrogen)/number_of_tanks; 
h_tank_first = volume_tank_first/(pi()*r_tank^2);
surface_area = number_of_tanks*(2*pi()*r_tank*h_tank_first)+(2*pi()*(r_tank^2));

%% Specs for tank not hydrogen
yield_stress = 4.61E8; %Pa
density_A = 2.8E3; %kg/m^3
SF = 2.5; 
stress_SF = SF*yield_stress;
storing_pressure_first = 101325;
thickness = ((storing_pressure_first)*(r_tank*2)*10)/(2*stress_SF); %10 for sanity
volume_tank = pi()*(h_tank_first+thickness)*(((r_tank+thickness)^2)-((r_tank)^2));
tank_weight = density_A*volume_tank;
end 
