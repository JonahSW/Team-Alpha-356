function [surface_area] = tank_size(total_mass_fraction, M_d, r_tank, number_of_tanks,str)
%Input

%Output

%% Calculating Tank Size (Cylinder)
M_w_first = total_mass_fraction*M_d; %Wet Mass [MT]
M_p_first = M_w_first - M_d;
M_P_first_kg = M_p_first *1000; 


if str == "hydrogen"
    liquid_hydrogen = 71; % kg/m^3
    volume_tank_first = (M_P_first_kg/liquid_hydrogen)/number_of_tanks;
    h_tank_first = volume_tank_first/(pi()*r_tank^2);
    surface_area_one = number_of_tanks*((2*pi()*r_tank*h_tank_first)+(2*pi()*(r_tank^2)));
    surface_area = surface_area_one;
elseif str == "xenon" 
    xenon_density = 3057; %kg/m^3
    volume_tank_first = (M_P_first_kg/xenon_density)/number_of_tanks; 
    radius_of_tank = (volume_tank_first*(3/4)*(1/pi))^(1/3);
    surface_area_one = 4*pi*radius_of_tank^2;
%     surface_area_one = 4*pi*(r_tank^2);
    surface_area = surface_area_one*number_of_tanks
else 
    disp("Can't Perform Calculation")
end 
    
    
%% Specs for tank not hydrogen
% yield_stress = 4.61E8; %Pa
% density_A = 2.8E3; %kg/m^3
% SF = 2.5; 
% stress_SF = SF*yield_stress;
% storing_pressure_first = 101325;
% thickness = ((storing_pressure_first)*(r_tank*2)*10)/(2*stress_SF); %10 for sanity
% volume_tank = pi()*(h_tank_first+thickness)*(((r_tank+thickness)^2)-((r_tank)^2));
% tank_weight = density_A*volume_tank;
end 
