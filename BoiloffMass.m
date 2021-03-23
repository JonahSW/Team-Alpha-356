function [BOM_mt, BOM_rate] = BoiloffMass(t_T,a_2,r_tank,number_of_tanks,total_mass_fraction_first,M_d)
%mdp81
%% Inputs
% t_T: Transfer Time (years)
% a_2: Distance to Ceres (AU)
% r_tank: Radius of Tank (m)
% number_of_tanks: Number of Tanks Used on Spacecraft (value)
% total_mass_faction: Mass fraction for half the mission (value)
% M_d = Dry Mass of Spacecraft (MT)
%% Outputs
% BOM_mt = Boil off Mass (MT)
% BOM_rate = rate of Boil off Mass per day (kg/second)
%% Calculating Temp as a function of Distance from the Sun 
M_w_first = total_mass_fraction_first*M_d; %Wet Mass [MT]
M_p_first = M_w_first - M_d;
a_u = 1.5E+8; % AU to Km
T_s = 5800; %temperature of the sun in K
R_s = 0.7E9; %Radius of the sun
distance_C = a_2 * a_u; % in km
distance_C_m = distance_C*1000;
distance_E_m = a_u *1000;
R_s_p = distance_E_m:100000000:distance_C_m;
Temp_2 = 2.7; %k
emissivity_MLI = (6.8E-4)*(Temp_2^(0.67)); % assumed to be the same for all number of layers
k_mli = (1E-5)*(Temp_2) ; %MLI conductivity 10^-5 W/mK

%Solar Constants
T_sun = 5772; %temperature of the sun in K
R_sun = 6.957e8; %Solar Equatorial Radius (m)
sigma = 5.67037e-8;% Stefan-Boltzmann Constant (W/m^2*K^4)
solar_surface_flux = sigma*T_sun^4; % Solar Radiation Flux at surface W/m2
solar_constant = 1.3608e3; % Solar Constant (W/m2) <- solar radiation flux at 1 AU

AU = 1.496e11;%AU in m
distance_from_sun = AU:(((a_2*AU)-AU)/2650):a_2*AU;

solar_flux = solar_surface_flux*(R_sun./distance_from_sun).^2;
%% Calculating Tank Size (Cylinder)
M_P_first_kg = M_p_first *1000; 
liquid_hydrogen = 71; % kg/m^3
viscosity_hydrogen = ((137E6)/10); % gram/cm sec 10^6 -> kg/m sec
C_p = 9.58*1000; %J/kg*k
kinematic_viscosity_hydrogen = viscosity_hydrogen/liquid_hydrogen;
% setting height of tank
% r_tank = 4.4; %meters
% number_of_tanks = 3;
volume_tank_first = (M_P_first_kg/liquid_hydrogen)/number_of_tanks; 
% volume_tank_second = (M_P_second_kg/liquid_hydrogen)/number_of_tanks; 
h_tank_first = volume_tank_first/(pi()*r_tank^2);
% h_tank_second = volume_tank_second/(pi()*r_tank^2);
surface_area = number_of_tanks*(2*pi()*r_tank*h_tank_first)+(2*pi()*(r_tank^2));
%% Calculating Thickness of Tank (Based on Aluminum 7075)
storing_pressure1 = []; 
BOM = [];
q_dprime_conduction = [];
q_dprime_convection = [];
q_dprime_radiation = [];
q_tot = [];
R = [];
F_g = [];
Ra = [];
Nu = [];
h = [];
molar_mass = [];
yield_stress = 4.61E8; %Pa
density_A = 2.8E3; %kg/m^3
SF = 2.5; 
alpha = 23.2E-6; % 1/C
stress_SF = SF*yield_stress;
storing_pressure = 101325; 
storing_pressure_first = 101325;
k_A = 190; % W/m
k_liquid_hydrogen = 102/1000; %W/mK
thermal_diffusivity = k_liquid_hydrogen/(liquid_hydrogen*C_p);
T_1_c = -252.778; % degrees Celsius
Temp_1 = T_1_c + 273.15; %T_1 in Kelvin 
number_of_layers = 20; % MLI half inch thick
thickness_of_MLI = 20*0.5*(1/39.37);
S_B_C = 5.670E-8; %stefan Boltzmann Constant W/m^2K^4
G = 6.67430E-11;
mass_earth = 5.972E24;
mass_ceres = 9.1E20;
latent_heat_evap_hydrogen = 461*1000; %KJ/kg
time_vector = 1:((t_T)/length(R_s_p))*3.154E7:t_T*3.154E7;
Temp_saturation = [];
thickness = ((storing_pressure_first)*(r_tank*2)*10)/(2*stress_SF); %10 for sanity
thickness_mm = thickness * 1000;
volume_tank = pi()*(h_tank_first+thickness)*(((r_tank+thickness)^2)-((r_tank)^2));
tank_weight = density_A*volume_tank;
pressure_gas_vap = [];
pressure_gas = [];
volume_gas = [];
Temp_boiloff = [];
%% Calculating Boiloff 
for n = 1:length(R_s_p)
    for k = 1:length(R_s_p)
        F_g(k) = (G* (M_w_first*1000))/(R_s_p(k)^2);
        Ra(k) = (alpha*(Temp_2-Temp_1)*(h_tank_first)*F_g(k))/(kinematic_viscosity_hydrogen*thermal_diffusivity);
        if Ra(k)<= 1E7
            Nu(k) = 0.642*((Ra(k))^(1/6));
        elseif Ra(k)<1E7 && Ra(k)<=1E10
            Nu(k) = 0.167*((Ra(k))^(1/4));
        else 
            Nu(k) = 0.00053*((Ra(k))^(1/2));
        end
        h(k) = (k_liquid_hydrogen/h_tank_first)*Nu(k);
        R(k) = (thickness/k_A)+(thickness_of_MLI/k_mli);
        q_dprime_conduction(k) = -(Temp_1-Temp_2)/R(k);
        q_dprime_convection(k) = -h(k)*(Temp_1-Temp_2); % not including space convection, dont know if that adds anything else
        q_dprime_radiation(k) = -emissivity_MLI*S_B_C*((Temp_1)^4 - (Temp_2)^4);
        q_tot(k) = surface_area*(q_dprime_conduction(k)+q_dprime_convection(k)+q_dprime_radiation(k)+solar_flux(k));
    end 
    Temp_saturation(n) = 1/((1/Temp_1) - ((287*log(storing_pressure/storing_pressure_first))/(latent_heat_evap_hydrogen)));
    pressure_gas_vap(n) = (287*(Temp_saturation(n))*q_tot(n)*time_vector(n))/(latent_heat_evap_hydrogen*volume_tank_first)-storing_pressure_first;
    pressure_gas(n) = (287*(Temp_saturation(n))*q_tot(n)*time_vector(n))/(latent_heat_evap_hydrogen*volume_tank_first);
    if storing_pressure>=storing_pressure_first
        storing_pressure = storing_pressure_first;
    else 
        storing_pressure = pressure_gas(n);
    end 
    volume_gas(n) = ((287*time_vector(n)*(Temp_saturation(n)-Temp_1))/pressure_gas_vap(n))-volume_tank_first;
    Temp_boiloff(n) = ((volume_gas(n)*pressure_gas(n))/C_p)-Temp_saturation(n);
    BOM(n) = (pressure_gas_vap(n)*volume_gas(n))/(287*Temp_boiloff(n));
end 
BOM_important = BOM(30:length(R_s_p));
time_start = time_vector(1);
new_time = [];
for j = 1:length(R)
    new_time(j) = time_start + time_vector(j);
end 
time_important = new_time(30:length(R_s_p));
BOM_Rate= real(BOM_important./time_important);
BOM_rate = mean(BOM_Rate);
Distance = distance_from_sun((30:length(R_s_p)));
R_S_P = (a_u*1000)*(R_s_p((30:length(R_s_p))));
limit = ((R_s_p(30))/(a_u*1000));
figure(2)
plot(R_S_P,BOM_Rate)
xlabel('R_{sp}')
ylabel('BOM rate')
xlim([limit t_T])
figure(3)
yyaxis left
plot(time_important/3.154E7,BOM_important)
ylabel('BOM')
hold on 
yyaxis right
plot(time_important/3.154E7,Distance)
ylabel('Distance (m)')
xlabel('Time (years)')
Boil_off_mass = real(sum(BOM_important)); %kg
BOM_mt = Boil_off_mass/1000; %MT
%% Accounting for the Boil Off Mass of Liquid Hydrogen
M_P_TOT = M_p_first + BOM_mt;
M_W_TOT = M_P_TOT+M_d;
M_P_per_tank = M_P_TOT/number_of_tanks;
BOM_mt_pertank = BOM_mt/number_of_tanks;
volume_tank_BOM = ((M_P_TOT*1000)/liquid_hydrogen)/number_of_tanks; 
h_tank_first_BOM = volume_tank_BOM/(pi()*r_tank^2);
end 
