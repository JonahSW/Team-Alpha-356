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
            target_v1 = delV_ppED(i)
            target_h1 = h_ED(i)
            target_mo_ma1 = mo_ma_ED(i)
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
            target_v2 = delV_ppCC(i)
            target_h2 = h_CC(i)
            target_mo_ma2 = mo_ma_CC(i)
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
            target_v3 = delV_ppCD(i)
            target_h3 = h_CD(i)
            target_mo_ma3 = mo_ma_CD(i)
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
            target_v4 = delV_ppEC(i)
            target_h4 = h_EC(i)
            target_mo_ma4 = mo_ma_EC(i)
            break
    end 
end 

%% total Del V 

target_delV = target_v1+target_v2+target_v3+target_v4
total_mass_fraction = target_mo_ma1*target_mo_ma2*target_mo_ma3*target_mo_ma4




%% Calculating Mass

M_d = 167.056; % Dry Mass [MT]
M_w = total_mass_fraction*M_d %Wet Mass [MT]
M_p = M_w - M_d % Propellant Mass [MT]
















