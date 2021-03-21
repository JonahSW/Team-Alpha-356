%mdp81, jjs280
%03/19/2021
%Code for estimating the boil off rate
%Takes a given heat flux and tank size and produces a boil off rate in kg/day
%Uses refprop to get hydrogen properties

function boil_off()

%% Solar System Constants
%Earth
earth_ephemeris = read_ephemeris(2);
a_earth = mean(earth_ephemeris(3,:));
in_earth = mean(earth_ephemeris(4,:));
earth_tilt = 23.43928*(pi/180);
launch_inclination = 28.573*(pi/180);
earth_radius = 6378000;
mu_earth = 3.986004418e14;% Earth Gravitational Parameter [m^3/sec^2] 
%Ceres
ceres_ephemeris = read_ephemeris(4);
a_ceres = mean(ceres_ephemeris(3,:));
in_ceres = mean(ceres_ephemeris(4,:));
ceres_tilt = 4*(pi/180);
ceres_radius = 469730;
mu_ceres = 6.26325e10;% Ceres Gravitational Parameter [m^3/sec^2]
%Sun
T_sun = 5772; %temperature of the sun in K
R_sun = 6.957e8; %Solar Equatorial Radius (m)

%% Calculating Temp as a function of Distance from the Sun 

distance_C = a_ceres * a_u; % in km
distance_C_m = distance_C*1000;
distance_E_m = a_u *1000;
R_s_p = distance_E_m:100000000:distance_C_m;

Temp_2 = [];
emissivity_MLI = [];
k_mli = [];
for k = 1:length(R_s_p)
    Temp_2(k) = T_sun*sqrt(R_sun/(2*R_s_p(k)));
    emissivity_MLI(k) = (6.8E-4)*(Temp_2(k)^(0.67)); % assumed to be the same for all number of layers
    k_mli(k) = (1E-5)*(Temp_2(k)) ; %MLI conductivity 10^-5 W/mK
end 

figure(1)
plot(R_s_p,Temp_2)
xlabel('R_{sp}')
ylabel('T_2')
title('Temperature as a function of distance away from the Sun')

%% Calculating Boiloff 
k_A = 190; % W/m
k_liquid_hydrogen = 102/1000; %W/mK
thermal_diffusivity = k_liquid_hydrogen/(liquid_hydrogen*C_p);
T_1_c = -252.778; % degrees Celsius
Temp_1 = T_1_c + 273.15; %T_1 in Kelvin 
number_of_layers = 20; % MLI half inch thick
thickness_of_MLI = 20*0.5*(1/39.37);
S_B_C = 5.670E-8; %stefan Boltzmann Constant W/m^2K^4
R = [];
q_dprime_conduction = [];
q_dprime_convection = [];
q_dprime_radiation = [];
q_tot = [];

mass_earth = 5.972E24;
mass_ceres = 9.1E20;

G = 6.67430E-11;
F_g = [];
Ra = [];
Nu = [];
h = [];
for k = 1:length(R_s_p)
    F_g(k) = (G* (M_w_first*1000))/(R_s_p(k)^2);
    Ra(k) = (alpha*(Temp_2(k)-Temp_1)*(h_tank_first)*F_g(k))/(kinematic_viscosity_hydrogen*thermal_diffusivity);
    if Ra(k)<= 1E7
        Nu(k) = 0.642*((Ra(k))^(1/6));
    elseif Ra(k)<1E7 && Ra(k)<=1E10
        Nu(k) = 0.167*((Ra(k))^(1/4));
    else 
        Nu(k) = 0.00053*((Ra(k))^(1/2));
    end
    h(k) = (k_liquid_hydrogen/h_tank_first)*Nu(k);
    R(k) = (thickness/k_A)+(thickness_of_MLI/k_mli(k));
    q_dprime_conduction(k) = -(Temp_1-Temp_2(k))/R(k);
    q_dprime_convection(k) = -h(k)*(Temp_1-Temp_2(k)); % not including space convection, dont know if that adds anything else
    q_dprime_radiation(k) = -emissivity_MLI(k)*S_B_C*((Temp_1)^4 - (Temp_2(k))^4);
    q_tot(k) = surface_area*(q_dprime_conduction(k)+q_dprime_convection(k)+q_dprime_radiation(k));
end 

figure(2)
plot(Ra,Nu)
xlabel('Ra')
ylabel('Nu')


latent_heat_evap_hydrogen = 461; %KJ/kg
BOM = [];
for k = 1:length(R_s_p)
    BOM(k) = q_tot(k)/latent_heat_evap_hydrogen; % kg
end 

figure(3)
plot(R_s_p,BOM)
xlabel('R_{sp}')
ylabel('BOM')

Boil_off_mass = sum(BOM) %kg
BOM_mt = Boil_off_mass/1000 %MT

%% Accounting for the Boil Off Mass of Liquid Hydrogen
M_P_TOT = M_p_first + BOM_mt
M_W_TOT = M_P_TOT+M_d

volume_tank_BOM = ((M_P_TOT*1000)/liquid_hydrogen)/number_of_tanks; 
h_tank_first_BOM = volume_tank_BOM/(pi()*r_tank^2)

end