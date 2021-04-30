% AAE 450 Fall 2001
% Vehicle propellant budget and tank sizing code
% Written by Casey Kirchner for Spring 2001, modified by Jon Edwards for Fall 2001

clear all;

%%%%%%%%%%%%%

% Constants %

%%%%%%%%%%%%%

g_0 = 9.81;  % m/s^2  This is a scaling factor to convert between force and mass, NOT necessarily to denote Earth's gravitational influence g_mars = 3.71;  % m/s^2  Acceleration of gravity on Mars

Ru = 8314.3;	% J/kmol-K, Universal Gas Constant

gamma = 1.66;  % Isentropic parameter for helium

Mw_He = 4.002602; % kg/kmol, molecular weight of helium

D_prop = 1.0;	% m, propellant tank diameter 
rho_g = 1550;  % kg/m^3 density of graphite (Table 5.15, Humble)

Ftu_g = 895000000;  % Pa, material strength of graphite (Table 5.15, Humble)

fs = 2;	% safety factor (Humble p. 269)

rho_al = 2800; % kg/m^3 (alloy 2219, Table 5.15, Humble)

rho_MMH = 0.878*1000;  % density of fuel - convert from g/cc (TEP output) to kg/m^3

rho_NTO = 1.444*1000;  % density of oxidizer

tw_al = 0.5/1000; % thickness of propellant tank aluminum liner  1/2 mm --> m


% Entry Angle Change Parameters
gamma_deg = 0; % deg, Angle Change
gamma = gamma_deg*(pi/180); % rad, change to radians
Ve = 8390; % m/s Speed at initial entry


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Vehicle Propulsion Calculations %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Engine Specifications

Isp_RCS = 306;% R-40A Marquardt RCS Engines 

Isp_retro = 313;% Shuttle OMS Engine (Hammond p.215)

landing_mass = 47.324;% tonnes Final mass as vehicle touches down on Mars surface

Ae_At = 55;  % Shuttle OMS expansion ratio

O_F = 1.6;  % O/F ratio Shuttle OMS Engine (Hammond p.215)


% Delta V required for spin-up maneuver
DV_tether_sul = 99.51;  % m/s (Table 4.4.2 Project PERForM Final Report)



% Delta V required for course corrections 

DV_RCS = 100*1.02;  % m/s (Spring 2001 Data, Table 4.4.2)


% Delta V required for entry angle change
DV_entry = Ve*sin(gamma); % m/s


% Delta V required for perigee lowering/raising in parking orbit 

DV_orbit = 0; % m/s   


% Apply rocket equation to solve for Hab mass fractions:  DV = g_0*Isp*ln(Mfinal/Minitial) - g_local*tb (neglect drag term);

% Parking orbit perigee burn phase

Mfraction_orbit = exp(-DV_orbit/(g_0*Isp_retro));

Minitial_orbit = landing_mass/Mfraction_orbit; % tonnes

prop_orbit_mass = Minitial_orbit - landing_mass; % tonnes

Mfinal_entry = Minitial_orbit; %tonnes
% Entry angle change burn phase
Mfraction_entry = exp(-DV_entry/(g_0*Isp_retro));
Minitial_entry = Mfinal_entry/Mfraction_entry; % tonnes
prop_entry_mass = Minitial_entry - Mfinal_entry; % tonnes
Mfinal_RCS = Minitial_entry; %tonnes

% Enroute maneuvering/attitude phase

Mfraction_RCS = exp(-DV_RCS/(g_0*Isp_RCS));

Minitial_RCS = Mfinal_RCS/Mfraction_RCS;

spun_mass = Mfinal_RCS;

prop_RCS_mass = Minitial_RCS - spun_mass;

Mfinal_sul = Minitial_RCS;

% Mass fractions for Spin-up phase
Mfraction_sul = exp(-DV_tether_sul/(g_0*Isp_RCS));

Minitial_sul = Mfinal_sul/Mfraction_sul;

prop_spin_up_mass = Minitial_sul - Mfinal_sul;

launched_mass = Minitial_sul;



%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Delta V and Prop Masses %

%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Total delta V required for Hab RCS/Retro

HAB_DV = DV_tether_sul + DV_RCS + DV_entry + DV_orbit; % m/s

% Total Hab propellant mass required

prop_total_mass = prop_orbit_mass + prop_entry_mass + prop_RCS_mass + prop_spin_up_mass; % metric tons


%%%%%%%%%%%%%%%%%%
% Engine Sizing  %
%%%%%%%%%%%%%%%%%%

% Thrust level

F_HAB = 26689;  % N, Thrust of Shuttle OMS engines

num_RCS = 24; % Number of RCS Marquardt R-40A engines

% Engine Physical Dimensions  
ME_HAB = 118; % kg (All OMS Engine Data from Astronautix.com http://www.astronautix.com/engines/ome.htm)
ME_RCS = 10 * num_RCS; % kg (All RCS Engine Data from Astronautix.com http://www.astronautix.com/engines/r40a.htm)

ME_HAB_total = ME_HAB + ME_RCS; % kg 

LE_HAB = 1.96; % m 

DE_HAB = 1.17;  % m 

LE_RCS = 1.0; % m
DE_RCS = 0.5; % m
total_engine_weight = ME_HAB_total; % kg

HAB_engine_weight = ME_HAB_total; % kg


%%%%%%%%%%%%%%%
% Tank Sizing %
%%%%%%%%%%%%%%%



% Propellant Masses

ox_total = (O_F/(O_F+1))*(prop_total_mass*1000);  % kg (Humble Eq 4.58)

fu_total = (1/(O_F+1))*(prop_total_mass*1000);  % kg (Humble Eq 4.58)

% Propellant Volumes

ox_volume = ox_total/rho_NTO*1.03; % m^3  Extra 3% for ullage (Humble p.268)

fu_volume = (4/3)*(3.45/2)^3*pi;%fu_total/rho_MMH*1.03; % m^3

% Propellant Tank Pressures

V_flow = 10;  % m/s typical flow velocity (Humble p.207)

Pc = 862000;  % Pa, Chamber pressure of OMS Engine (from http://www.astronautix.com/engines/ome.htm)

DP_feed = 50000;  % Pa, Max pressure drop through feedsystem (Humble Eq 5.18)
DP_cool = 0.15*Pc; % Pa, Pressure loss through cooling system (Humble Eq 5.19)

DP_dynamic_ox = 0.5*rho_NTO*(V_flow)^2;  % Oxidizer dynamic pressure drop (Humble Eq 5.16)

DP_dynamic_fu = 0.5*rho_MMH*(V_flow)^2; % Fuel dynamic pressure drop (Humble Eq 5.16)

DP_injector = 0.3*Pc;  % Delta P through injector for throttled system (Humble Eq 5.21)

DP_total = DP_feed + DP_dynamic_ox + DP_dynamic_fu + DP_injector + DP_cool; % Pa, Total Pressure loss through feed system

Ptank = 8.62e6;%Pc + DP_total; % Pa, Tank Pressure needed for OMS Chamber Pressure to be reached (Humble Eq 5.24)

% Pressurant tank pressures

Ti = 300;  % K, temp of pressurant tank before blowdown (Humble p.278)

Tf = 200;  % K, temp of pressurant tank after blowdown (Humble p.278)

HAB_press_vol = ox_volume + fu_volume; % Initial guess of final volume of pressurant required (Humble p. 279)

HAB_press_mass = Ptank*HAB_press_vol*Mw_He/(Ru*Tf); % kg, mass of the pressurant required (Humble Eq 5.83)

press_on_HAB_vol = HAB_press_mass*Ru*Ti/(Ptank*Mw_He); % m^3, Required pressurant tank volume (Ideal Gas Law)

burst_pressure = fs * Ptank; % Structure's design burst pressure (Humble Eq 5.73)
% Tank Lengths

Loxt_HAB = (ox_volume - (4/3)*pi*(D_prop/2)^3)/(pi*(D_prop/2)^2); % m, Cylindrical length of ox tank (Humble Eq 5.78)

Lfut_HAB = (fu_volume - (4/3)*pi*(D_prop/2)^3)/(pi*(D_prop/2)^2); % m, Cylindrical lenght of fuel tank (Humble Eq 5.78)

Dprt_HAB = 2*(3*press_on_HAB_vol/(4*pi))^(1/3);  % m, Diameter of Pressurant Tank (Humble Eq 5.74)

% Tank Surface Areas

Soxt_HAB = pi*D_prop^2 + pi*D_prop*Loxt_HAB; % m^2, Surface area of ox tank (Humble Eq 5.79 and Eq 5.75)

Sfut_HAB = pi*D_prop^2 + pi*D_prop*Lfut_HAB; % m^2, Surface area of fuel tank (Humble Eq 5.79 and Eq 5.75)

Sprt_HAB = pi*Dprt_HAB^2; % m^2, Surface area of pressurant tank (Humble Eq 5.75)

% Tank Wall Thicknesses

tw_oxt_HAB = (D_prop*burst_pressure)/(2*Ftu_g); % m, ox tank thickness (Humble Eq 5.80)

tw_fut_HAB = (D_prop*burst_pressure)/(2*Ftu_g); % m, fuel tank thickness (Humble Eq 5.80)

tw_prt_HAB = (Dprt_HAB*burst_pressure)/(4*Ftu_g); % m, pressurant tank thickness (Humble Eq 5.76)

% Tank Masses

Moxt_HAB = rho_g*Soxt_HAB*tw_oxt_HAB + rho_al*Soxt_HAB*tw_al; % kg, ox tank mass (includes aluminum liner)

Mfut_HAB = rho_g*Sfut_HAB*tw_fut_HAB + rho_al*Sfut_HAB*tw_al; % kg, fuel tank mass (includes aluminum liner)

Mprt_HAB = rho_g*Sprt_HAB*tw_prt_HAB; % kg, pressurant tank mass

total_liner_mass = rho_al*tw_al*(Soxt_HAB + Sfut_HAB);  % kg, total mass of tank liners

HAB_tank_mass = Moxt_HAB + Mfut_HAB + Mprt_HAB;   % kg, total mass of tanks



%%%%%%%%%%%%%%

% Inert mass %

%%%%%%%%%%%%%%



% According to the AAE 439/539 text, we add 10% of the total inert mass for support structure (p. 281 Humble)

HAB_support_mass = 0.10*(ME_HAB_total + HAB_tank_mass);

HAB_inert_mass = ME_HAB_total + HAB_tank_mass + HAB_support_mass;



%%%%%%%%%%

% Output %

%%%%%%%%%%



fprintf(1,'%1s\n',' ');

fprintf(1,'%19s\n','Delta V Budget, m/s');

fprintf(1,'%1s\n',' ');

fprintf(1,'%32s %8.2f\n','Mars Orbit DV:',DV_orbit);
fprintf(1,'%32s %8.2f\n','Entry Angle Change DV:',DV_entry);

fprintf(1,'%32s %8.2f\n','Spin-up DV',DV_tether_sul);

fprintf(1,'%32s %8.2f\n','Enroute RCS/Maneuvering DV:',DV_RCS);

fprintf(1,'%32s %8.2f\n','Total DV:',HAB_DV);

fprintf(1,'%1s\n',' ');

fprintf(1,'%6s\n','Masses');

fprintf(1,'%32s %15s\n',' ','  Hab Subsystem');

fprintf(1,'%32s %15.2f\n','OMS Engine mass (kg):       ',ME_HAB);

fprintf(1,'%32s %15.2f\n','RCS Engine mass (total) (kg):    ',ME_RCS);

fprintf(1,'%32s %15.2f\n','Oxidizer tank masses (kg):      ',Moxt_HAB);

fprintf(1,'%32s %15.2f\n','Fuel tank masses (kg):          ',Mfut_HAB);

fprintf(1,'%32s %15.2f\n','Pressurant tank masses (kg):    ',Mprt_HAB);

fprintf(1,'%32s %15.2f\n','Struct. support mass (kg):      ',HAB_support_mass);
fprintf(1,'%1s\n',' ');

fprintf(1,'%32s %15.2f\n','Total inert mass (kg):          ',HAB_inert_mass);

fprintf(1,'%1s\n',' ');

fprintf(1,'%32s %15.2f\n','Oxidizer masses (kg):           ',ox_total);

fprintf(1,'%32s %15.2f\n','Fuel masses (kg):               ',fu_total);

fprintf(1,'%32s %15.2f\n','Pressurant masses (kg):         ',HAB_press_mass);

fprintf(1,'%32s %15.2f\n','Total prop/press masses (kg):   ',ox_total+fu_total+HAB_press_mass);

fprintf(1,'%1s\n',' ');

fprintf(1,'%32s %15.2f\n','Total mass (kg):        ',ox_total+fu_total+HAB_press_mass+HAB_inert_mass);

fprintf(1,'%1s\n',' ');

fprintf(1,'%15s\n','Tank and Engine Geometry');

fprintf(1,'%32s %15s %15s\n','                                ','  Length (m)','  Diameter (m)');

fprintf(1,'%32s %15.2f %15.2f\n','OMS Engine:                   ',LE_HAB,DE_HAB);

fprintf(1,'%32s %15s %15.2f\n','OMS Engine Throat:            ','            N/A',(pi*DE_HAB^2/4)/Ae_At);

fprintf(1,'%32s %15.2f %15.2f\n','RCS Engines on Hab :',LE_RCS,DE_RCS);

fprintf(1,'%32s %15.2f %15.2f\n','Oxidizer tank on Hab:           ',Loxt_HAB+D_prop,D_prop);

fprintf(1,'%32s %15.2f %15.2f\n','Fuel tank on Hab:               ',Lfut_HAB+D_prop,D_prop);

fprintf(1,'%32s %15s %15.2f\n','Pressurant tank on Hab:         ','      Spherical',Dprt_HAB);

fprintf(1,'%1s\n',' ');



