function [distance_from_sun,Boil_off_rate,BOM,BOM_duration] = BoiloffMasslockheed(str,a_1,a_2,t_T,surface_area)

%% Input

%% Constants
a_earth = 1;
a_ceres = 2.766; 
a_u = 1.60E11; %AU to m
spacing = 1000;
years_to_secs = 3.15E7;
C_c = 8.85E-5;
C_r = 5.39E-7;
c_p_xenon = 0.158 *1000; %J/kg*K
T_crit_xenon = 289.65; % K
n = 3; 
N = 20; % 20 layers/cm
percent = 0.4; %Percent of the Spacecraft seeing the solar flux
emissivity_MLI = 0.0331;
%% Vectors for distance and time
if a_1 == a_2
    distance_from_sun = a_1*a_u*ones(1,spacing);
else 
    distance_from_sun = a_1*a_u:(((a_2-a_1)*a_u)/spacing):a_2*a_u;
end 
time_vector = ((t_T*years_to_secs)/length(distance_from_sun))*ones(1,length(distance_from_sun));
%%
if str == "hydrogen"
    T_1 = 20.372; % hydrogen (k)
    L_vap = 461*1000; %J/kg
    [T_2] = EffectiveTemp(a_1,a_2,a_earth,a_ceres,distance_from_sun,surface_area,T_1);
    T_m = (T_2+T_1)./2;
    q = ((((C_c * N^(2.56) * T_m).*(T_2-T_1))./n) + ((C_r * emissivity_MLI).*((T_2.^(4.67))-(T_1^(4.67)))./n))./1000;
elseif str == "xenon"
    T_1 = 120; % xenon (k)
    [T_2] = EffectiveTemp(a_1,a_2,a_earth,a_ceres,distance_from_sun,surface_area,T_1);
    L_vap = (96.29*1000)+(c_p_xenon.*(T_crit_xenon - T_1));
    T_m = (T_2+T_1)./2;
    q = ((((C_c * N^(2.56) * T_m).*(T_2-T_1))./n) + ((C_r * emissivity_MLI).*((T_2.^(4.67))-(T_1^(4.67)))./n))./1000;
else 
    disp("Can't Perform Calculation") 
end 
flux = (percent*q.*surface_area);
Boil_off_rate = flux./L_vap;
BOM = Boil_off_rate.*time_vector;
BOM_duration = (sum(BOM))/1000;
end

   
    