clear 
close all
clc
%% Constants 
a_1 = 1; % Distance from Earth to Sun 
a_2 = 2.766; % Distance from Ceres to Sun 
u_O = 1.33E+11; % Gravitation Parameter Sun (km^3/s^2)
a_u = 1.5E+8; % AU to Km
mu_A = 398600; % Gravitational Parameter of Earth (km^3/s^2)
mu_B = 62.630; % Gravitational Parameter of Ceres (km^3/s^2)
r_OE = 6378; % Radius of Earth Km
r_OC = 469.730; % Radius of Ceres Km
EAI = 23.450; % Earth Axis of Inclination (Degrees)
COP = 10.110; % Ceres Orbit Plane (Degrees)
EPOI = 28.573; % Earth Parking Orbit Inclination (Degrees)
KSC = 27.573; % Kennedy Space Center (Degrees)
CAI = 4.00; % Ceres Axis Inclination (Degrees) 
g = 9.807; % gravity m/s^2
DPC = (EAI +COP) - EPOI;
constants = ["Orbital Height"; "Radius + Orbital Height";"Circular Orbit Velocity";"Ecentricity of Ellipse";"Semi Minor Axis";"Flight Path Angle at Second Burn";"Departure/Capture Velo";"Velo Necessary to Match Circular Velo";"Velocity Associated with Plane Change";"DelV then Plan Change";"DelV and Plane Change at the Same Time";"Mass Fraction"];
%% Input Parameter
h_ED = 200:1:500; %Input Height (Earth Departure). Iterating to see multiple values
h_CC = 100:1:400; %Input Height (Ceres Capture). Iterating to see multiple values
% h_ED = 408;
% h_CC = 100;
h_EC = h_ED; %Input Height (Earth Capture). Iterating to see multiple values
h_CD = h_CC; %Input Height (Ceres Departure). Iterating to see multiple values
l_1 = length(h_ED); 
l_2 = length (h_CC);
%% Earth-Ceres Hohmann Transfer
V_A = sqrt(u_O/a_u); % Va in km/s
V_B = sqrt(u_O/(a_u * a_2)); % Vb in km/s
delV_A = V_A * (((2*a_2/(a_1+a_2))^(1/2))-1); % First Burn in km/s 
delV_B = V_B * (1-((2*a_1/(a_1+a_2))^(1/2))); % First Burn in km/s 
delV_tot = delV_A + delV_B; % total burn velocity in km/s 
e_T = (((delV_A/V_A)+1)^2)-1; % Ecentricity of Transfer
a_T = (a_2+a_1)/2; % Semi Major Axis
t_T = 0.5*sqrt(a_T^3); % transfer time (years)
T_1 = 1; % Orbit Period of the Earth (years)
T_2 = sqrt((a_2/a_1)^3); % Orbit Period of Ceres
S = 1/((1/T_1)-(1/T_2)); % Synodic Period
alpha_pi = (((a_1+a_2)/(2*a_1))^(3/2))-1; % time used for wait time in years
y_w = 1:1:5; %iterating years to wait
for i = 1:5
    t_w = (y_w(i) - alpha_pi)*S; % wait time on ceres (years)
    if t_w <0
        i = i+1;
    else
        break
    end   
end 
t_t = (2*t_T) + t_w; %Total Mission time (years)
%% Earth Departure- Parameters and Solve
a_OED = []; % Radius + Height (km)
V_cirED = []; % Circular Orbit Velo (km/s)
a_ED = mu_A/(delV_A^2); % Semi Major Axis of Ellipse 
e_ED = []; % Ecentricity of Ellipse
b_ED = []; % Semi Minor Axis
psi_ED = []; %flight path angle at 2nd burn (degrees)
V_HED = []; % Departure Velo (km/s)
delV_DED = []; % Velocity necessary to match Circular Velo (km/s)
delV_pcED = []; % Velo Necessary for Plane Change (km/s)
delV_pED = []; % Del V then plane change (km/s)
delV_ppED =[]; % del V and plane change at the same time (km/s)
ISP = 950; % Input ISP Value (s)
mo_ma_ED = []; % mass fraction
for i = 1:l_1
    a_OED(i) = h_ED(i) + r_OE;
    V_cirED(i) = sqrt(mu_A/a_OED(i));
    e_ED(i) = (a_OED(i)/a_ED)+1; 
    b_ED(i) = a_ED*sqrt(((e_ED(i))^2)-1);
    psi_ED(i) = atand(sqrt(((e_ED(i))^2)-1));
    V_HED(i) = (((e_ED(i)+1)/(e_ED(i)-1))^0.5)*delV_A;
    delV_DED(i) = V_HED(i) - V_cirED(i);
    delV_pcED(i) = 2*V_cirED(i)*sind(DPC/2);
    delV_pED(i) = delV_pcED(i) + delV_DED(i);
    delV_ppED(i) = sqrt(((V_cirED(i))^2)+((V_HED(i))^2) - (2*V_cirED(i)*V_HED(i)*cosd(DPC)));
    mo_ma_ED(i) = exp((1000*delV_DED(i))/(g*ISP));     
end 
matrix_ED = [constants [h_ED;a_OED;V_cirED;e_ED;b_ED;psi_ED;V_HED;delV_DED;delV_pcED;delV_pED;delV_ppED;mo_ma_ED]];
for i = 1:l_1
    if min(mo_ma_ED)==mo_ma_ED(i)
            target_v1 = delV_ppED(i);
            target_h1 = h_ED(i);
            target_mo_ma1 = mo_ma_ED(i);
            break
    end 
end      
%% Ceres Capture- Parameters and Solve
a_OCC = []; % Radius + Height (km)
V_cirCC = []; % Circular Orbit Velo (km/s)
a_CC = mu_B/(delV_B^2); % Semi Major Axis of Ellipse 
e_CC = []; % Ecentricity of Ellipse
b_CC = []; % Semi Minor Axis
psi_CC = []; %flight path angle at 2nd burn (degrees)
V_HCC = []; % Departure Velo (km/s)
delV_CCC = []; % Velocity necessary to match Circular Velo (km/s)
delV_pcCC = []; % Velo Necessary for Plane Change (km/s)
delV_pCC = []; % Del V then plane change (km/s)
delV_ppCC =[]; % del V and plane change at the same time (km/s)
mo_ma_CC = []; % mass fraction
check1 = [];
for i = 1:l_2
    a_OCC(i) = h_CC(i) + r_OC;
    V_cirCC(i) = sqrt(mu_B/a_OCC(i));
    e_CC(i) = (a_OCC(i)/a_CC)+1; 
    b_CC(i) = a_CC*sqrt(((e_CC(i))^2)-1);
    psi_CC(i) = atand(sqrt(((e_CC(i))^2)-1));
    V_HCC(i) = (((e_CC(i)+1)/(e_CC(i)-1))^0.5)*delV_B;
    delV_DCC(i) = V_HCC(i) - V_cirCC(i);
    delV_pcCC(i) = 2*V_cirCC(i)*sind(CAI/2);
    delV_pCC(i) = delV_pcCC(i) + delV_DCC(i);
    delV_ppCC(i) = sqrt(((V_cirCC(i))^2)+((V_HCC(i))^2) - (2*V_cirCC(i)*V_HCC(i)*cosd(CAI)));
    mo_ma_CC(i) = exp((1000*delV_DCC(i))/(g*ISP));  
end 
matrix_CC = [constants [h_CC;a_OCC;V_cirCC;e_CC;b_CC;psi_CC;V_HCC;delV_DCC;delV_pcCC;delV_pCC;delV_ppCC;mo_ma_CC]];
for i = 1:l_2
    if min(mo_ma_CC)==mo_ma_CC(i)
            target_v2 = delV_ppCC(i);
            target_h2 = h_CC(i);
            target_mo_ma2 = mo_ma_CC(i);
            break
    end 
end 
%% Ceres-Earth Hohmann Transfer
V_A2 = sqrt(u_O/(a_u*a_2)); % Va in km/s
V_B2 = sqrt(u_O/(a_u)); % Vb in km/s
delV_A2 = V_A2 * (((2*a_2/(a_1+a_2))^(1/2))-1); % First Burn in km/s 
delV_B2 = V_B2 * (1-((2*a_1/(a_1+a_2))^(1/2))); % First Burn in km/s 
delV_tot2 = delV_A2 + delV_B2; % total burn velocity in km/s 
e_T2 = (((delV_A2/V_A2)+1)^2)-1; % Ecentricity of Transfer
a_T2 = (a_2+a_1)/2; % Semi Major Axis
t_T2 = 0.5*sqrt(a_T2^3); % transfer time (years)
T_12 = 1; % Orbit Period of the Earth (years)
T_22 = sqrt((a_2/a_1)^3); % Orbit Period of Ceres
S2 = 1/((1/T_12)-(1/T_22)); % Synodic Period
alpha_pi2 = (((a_1+a_2)/(2*a_1))^(3/2))-1; % time used for wait time in years
y_w2 = 1:1:5; %iterating years to wait
for i = 1:5
    t_w2 = (y_w(i) - alpha_pi)*S; % wait time on ceres (years)
    if t_w2 <0
        i = i+1;
    else
        break
    end 
end 
t_t2 = (2*t_T2) + t_w2; %Total Mission time (years)
%% Ceres Departure 
a_OCD = []; % Radius + Height (km)
V_cirCD = []; % Circular Orbit Velo (km/s)
a_CD = mu_B/(delV_A2^2); % Semi Major Axis of Ellipse 
e_CD = []; % Ecentricity of Ellipse
b_CD = []; % Semi Minor Axis
psi_CD = []; %flight path angle at 2nd burn (degrees)
V_HCD = []; % Departure Velo (km/s)
delV_CCD = []; % Velocity necessary to match Circular Velo (km/s)
delV_pcCD = []; % Velo Necessary for Plane Change (km/s)
delV_pCD = []; % Del V then plane change (km/s)
delV_ppCD =[]; % del V and plane change at the same time (km/s)
mo_ma_CD = []; % mass fraction
for i = 1:l_2
    a_OCD(i) = h_CD(i) + r_OC;
    V_cirCD(i) = sqrt(mu_B/a_OCD(i));
    e_CD(i) = (a_OCD(i)/a_CD)+1; 
    b_CD(i) = a_CD*sqrt(((e_CD(i))^2)-1);
    psi_CD(i) = atand(sqrt(((e_CD(i))^2)-1));
    V_HCD(i) = (((e_CD(i)+1)/(e_CD(i)-1))^0.5)*delV_A2;
    delV_DCD(i) = V_HCD(i) - V_cirCD(i);
    delV_pcCD(i) = 2*V_cirCD(i)*sind(CAI/2);
    delV_pCD(i) = delV_pcCD(i) + delV_DCD(i);
    delV_ppCD(i) = sqrt(((V_cirCD(i))^2)+((V_HCD(i))^2) - (2*V_cirCD(i)*V_HCD(i)*cosd(CAI)));
    mo_ma_CD(i) = exp((1000*delV_DCD(i))/(g*ISP));    
end 
matrix_CD = [constants [h_CD;a_OCD;V_cirCD;e_CD;b_CD;psi_CD;V_HCD;delV_DCD;delV_pcCD;delV_pCD;delV_ppCD;mo_ma_CD]];
for i = 1:l_1
    if min(mo_ma_CD)==mo_ma_CD(i)
            target_v3 = delV_ppCD(i);
            target_h3 = h_CD(i);
            target_mo_ma3 = mo_ma_CD(i);
            break
    end 
end 
%% Earth Capture- Parameters and Solve
a_OEC = []; % Radius + Height (km)
V_cirEC = []; % Circular Orbit Velo (km/s)
a_EC = mu_A/(delV_B2^2); % Semi Major Axis of Ellipse 
e_EC = []; % Ecentricity of Ellipse
b_EC = []; % Semi Minor Axis
psi_EC = []; %flight path angle at 2nd burn (degrees)
V_HEC = []; % Departure Velo (km/s)
delV_DEC = []; % Velocity necessary to match Circular Velo (km/s)
delV_pcEC = []; % Velo Necessary for Plane Change (km/s)
delV_pEC = []; % Del V then plane change (km/s)
delV_ppEC =[]; % del V and plane change at the same time (km/s)
mo_ma_EC = []; % mass fraction
for i = 1:l_1
    a_OEC(i) = h_EC(i) + r_OE;
    V_cirEC(i) = sqrt(mu_A/a_OEC(i));
    e_EC(i) = (a_OEC(i)/a_EC)+1; 
    b_EC(i) = a_EC*sqrt(((e_EC(i))^2)-1);
    psi_EC(i) = atand(sqrt(((e_EC(i))^2)-1));
    V_HEC(i) = (((e_EC(i)+1)/(e_EC(i)-1))^0.5)*delV_B2;
    delV_DEC(i) = V_HEC(i) - V_cirEC(i);
    delV_pcEC(i) = 2*V_cirEC(i)*sind(DPC/2);
    delV_pEC(i) = delV_pcEC(i) + delV_DEC(i);
    delV_ppEC(i) = sqrt(((V_cirEC(i))^2)+((V_HEC(i))^2) - (2*V_cirEC(i)*V_HEC(i)*cosd(DPC)));
    mo_ma_EC(i) = exp((1000*delV_DEC(i))/(g*ISP));
end 
matrix_EC = [constants [h_EC;a_OEC;V_cirEC;e_EC;b_EC;psi_EC;V_HEC;delV_DEC;delV_pcEC;delV_pEC;delV_ppEC;mo_ma_EC]];
for i = 1:l_1
    if min(mo_ma_EC)==mo_ma_EC(i)
            target_v4 = delV_ppEC(i);
            target_h4 = h_EC(i);
            target_mo_ma4 = mo_ma_EC(i);
            break
    end 
end 
%% total Del V 
target_delV = target_v1+target_v2+target_v3+target_v4;
total_mass_fraction_first = target_mo_ma1*target_mo_ma2;
total_mass_fraction_second = target_mo_ma3*target_mo_ma4;
%% Calculating Mass
M_d = 167.056; % Dry Mass [MT]
M_w_first = total_mass_fraction_first*M_d; %Wet Mass [MT]
M_w_second = total_mass_fraction_first*M_d;
M_w = M_w_first+M_w_second;
M_p_first = M_w_first - M_d; % Propellant Mass [MT]
M_p_second = M_w_second - M_d;
M_P_first_kg = M_p_first*1000;
M_P_second_kg = M_p_second*1000;
%% Calculating Temp as a function of Distance from the Sun 

T_s = 5800; %temperature of the sun in K
R_s = 0.7E9; %Radius of the sun
distance_C = a_2 * a_u; % in km
distance_C_m = distance_C*1000;
distance_E_m = a_u *1000;
R_s_p = distance_E_m:100000000:distance_C_m;
Temp_2 = [];
emissivity_MLI = [];
k_mli = [];
for k = 1:length(R_s_p)
    Temp_2(k) = T_s*sqrt(R_s/(2*R_s_p(k)));
    emissivity_MLI(k) = (6.8E-4)*(Temp_2(k)^(0.67)); % assumed to be the same for all number of layers
    k_mli(k) = (1E-5)*(Temp_2(k)) ; %MLI conductivity 10^-5 W/mK
end 

figure(1)
plot(R_s_p,Temp_2)
xlabel('R_{sp}')
ylabel('T_2')
title('Temperature as a function of distance away from the Sun')
%% Calculating Tank Size (Cylinder)
liquid_hydrogen = 71; % kg/m^3
viscosity_hydrogen = ((137E6)/10); % gram/cm sec 10^6 -> kg/m sec
C_p = 9.58*1000; %J/kg*k
kinematic_viscosity_hydrogen = viscosity_hydrogen/liquid_hydrogen;
% setting height of tank
r_tank = 4.5; %meters
number_of_tanks = 3;
volume_tank_first = (M_P_first_kg/liquid_hydrogen)/number_of_tanks; 
volume_tank_second = (M_P_second_kg/liquid_hydrogen)/number_of_tanks; 
h_tank_first = volume_tank_first/(pi()*r_tank^2);
h_tank_second = volume_tank_second/(pi()*r_tank^2);

surface_area = 3*(2*pi()*r_tank*h_tank_first)+(2*pi()*(r_tank^2));

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
pressure_gas = [];

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
latent_heat_evap_hydrogen = 461; %KJ/kg
%% Calculating Boiloff 
for n = 1:length(R_s_p)
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

    for k = 1:length(R_s_p)
        BOM(k) = q_tot(k)/latent_heat_evap_hydrogen; % kg
    end 

    pressure_gas(n) = (BOM(n)* 287 * Temp_1)/volume_tank;
    storing_pressure1 = [storing_pressure pressure_gas];

    storing_pressure = storing_pressure + pressure_gas(n);
    thickness = ((storing_pressure)*(r_tank*2)*10)/(2*stress_SF); %10 for sanity
    thickness_mm = thickness * 1000;
    volume_tank = pi()*(h_tank_first+thickness)*(((r_tank+thickness)^2)-((r_tank)^2));
    tank_weight = density_A*volume_tank;
end 

figure(2)
plot(Ra,Nu)
xlabel('Ra')
ylabel('Nu')

figure(3)
plot(R_s_p,BOM)
xlabel('R_{sp}')
ylabel('BOM')


Boil_off_mass = sum(BOM); %kg


BOM_mt = Boil_off_mass/1000; %MT

%% Accounting for the Boil Off Mass of Liquid Hydrogen
M_P_TOT = M_p_first + BOM_mt
M_W_TOT = M_P_TOT+M_d

volume_tank_BOM = ((M_P_TOT*1000)/liquid_hydrogen)/number_of_tanks; 
h_tank_first_BOM = volume_tank_BOM/(pi()*r_tank^2)


