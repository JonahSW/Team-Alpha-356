%% Clean State
clear; close all; clc;
%% Constant
a_Earth = 1; % Distance from Earth to Sun 
a_Ceres = 2.766; % Distance from Ceres to Sun
r_Earth = 6378; % Radius of Earth (km)
r_Ceres = 469.730; % Radius of Ceres (km)
mu_Earth = 398600; % Gravitational Parameter of Earth (km^3/s^2)
mu_Ceres = 62.630; % Gravitational Parameter of Ceres (km^3/s^2)
Ba_Earth = 0.306; % Bond albedo of Earth
Ba_Ceres = 0.4; % Bond albedo of Ceres
T_Sun = 5772; % Temperture of the Surface of the Sun (K);
R_Sun = 696340; % Radius of the Sun (km);
ISP = 950;

%% Transfer
% From Earth to Ceres
% Hohamman Tranfer
[t_trans_EC,delV_tot_EC,delV_A_EC,delV_B_EC] = Hohmann(a_Earth, a_Ceres);
% Earth Departure
h_ED = 200:1:500; %Input Height
[MassFraction_ED,Target_v_ED,Target_h_ED,Design_ED] = DepCap(r_Earth,h_ED,mu_Earth,delV_A_EC,ISP);
% Ceres Capture
h_CC = 100:1:400; %Input Height
[MassFraction_CC,Target_v_CC,Target_h_CC,Design_CC] = DepCap(r_Ceres,h_CC,mu_Ceres,delV_B_EC,ISP);

% Ceres to Earth Wait
[t_wait] = Wait(a_Earth, a_Ceres); % Waiting Time on Ceres

% From Ceres to Earth
% Hohamman Tranfer
[t_trans_CE,delV_tot_CE,delV_A_CE,delV_B_CE] = Hohmann(a_Ceres, a_Earth);
% Ceres Departure
h_CD = h_CC; %Input Height
[MassFraction_CD,Target_v_CD,Target_h_CD,Design_CD] = DepCap(r_Ceres,h_CD,mu_Ceres,delV_A_CE,ISP);
% Earth Capture
h_EC = h_ED; %Input Height
[MassFraction_EC,Target_v_EC,Target_h_EC,Design_EC] = DepCap(r_Earth,h_EC,mu_Earth,delV_B_CE,ISP);

% Transfer Sumup
Target_h_E = Target_h_ED; % Orbit Hight of Earth
Target_h_C = Target_h_CC; % Orbit High of Ceres
Target_v_E = Target_v_ED; % Orbit Velocity of Earth
Target_v_C = Target_v_CC; % Orbit Velocity of Ceres
Target_delV = Target_v_ED + Target_v_CC + Target_v_CD + Target_v_EC;
% MassFraction_ED*MassFraction_CC = MassFraction_CD*MassFraction_EC;
MassFraction = MassFraction_ED*MassFraction_CC;

%% Propellant Mass Calculation
M_d = 167.056; % Dry Mass (MT)
M_w = MassFraction*M_d; % Wet Mass (MT)
M_ph = M_w - M_d; % Half-transfer Propellant Mass(MT)
M_p = M_ph*2;

%% Propellant Tank Calculation
% Asumptions on Propellant Tanks
r_out = 1.4; % Out Radius of Tanks
n = 3; % Number of Tanks
[Tank_A_out] = TankDimensions(r_out,n,M_p); % Outter Surface Area of All Tanks


%% Boiloff Rate





