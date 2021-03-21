%mdp81, jjs280
%03/21/2021
%Code for estimating the tank size and mass for cryogenic propellant storage
%Given a target wet mass, estimates the size and wall thickness of the tanks
%Tanks modelled as a cylinder with hemispherical endcaps
%Uses refprop to get hydrogen properties

function tank_size()
%% Calculating Tank Size (Cylinder with hemispherical endcaps)
liquid_hydrogen = 71; % kg/m^3
viscosity_hydrogen = ((137E6)/10); % gram/cm sec 10^6 -> kg/m sec
C_p = 9.58*1000; %J/kg*k
kinematic_viscosity_hydrogen = viscosity_hydrogen/liquid_hydrogen;
% setting height of tank
r_tank = 5; %meters
number_of_tanks = 3;
volume_tank_first = (M_P_first_kg/liquid_hydrogen)/number_of_tanks; 
volume_tank_second = (M_P_second_kg/liquid_hydrogen)/number_of_tanks; 
h_tank_first = volume_tank_first/(pi()*r_tank^2);
h_tank_second = volume_tank_second/(pi()*r_tank^2);

surface_area = 3*(2*pi()*r_tank*h_tank_first)+(2*pi()*(r_tank^2));

%% Calculating Thickness of Tank (Based on Aluminum 7075)
yield_stress = 4.61E8; %Pa
density_A = 2.8E3; %kg/m^3
SF = 2.5; 
alpha = 23.2E-6; % 1/C
stress_SF = SF*yield_stress;
storing_pressure = 101325; 
thickness = ((storing_pressure)*(r_tank*2)*10)/(2*stress_SF); %10 for sanity
thickness_mm = thickness * 1000;

volume_tank = pi()*(h_tank_first+thickness)*(((r_tank+thickness)^2)-((r_tank)^2));
tank_weight = density_A*volume_tank;


end