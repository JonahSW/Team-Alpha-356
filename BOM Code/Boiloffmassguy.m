function [BOM_important, time_important, BOM_mt] = Boiloffmassguy(t_T, R_s_p, q_tot,volume_tank_first, R, number_of_tanks)
%Input

%Output

T_1_c = -252.778; % degrees Celsius
Temp_1 = T_1_c + 273.15; %T_1 in Kelvin 
latent_heat_evap_hydrogen = 461*1000; %J/kg
t_T_sec = t_T*3.154E7;
time_vector = 0:((t_T_sec)/length(R_s_p)):t_T_sec;
R_constant = 4124.2; 
storing_pressure = 101325;
storing_pressure_first = 101325;
C_p = 9.58*1000; %J/kg*k
liquid_hydrogen = 71; % kg/m^3
BOM = [];
Temp_saturation = [];
pressure_gas_vap = [];
pressure_gas = [];
volume_gas = [];
Temp_boiloff = [];
for n = 1:length(R_s_p)
    Q_tot(n) = q_tot(n);
    Temp_saturation(n) = 1/((1/Temp_1) - ((R_constant*log(storing_pressure/storing_pressure_first))/(latent_heat_evap_hydrogen)));
    pressure_gas(n) = (R_constant*Temp_saturation(n)*Q_tot(n)*(time_vector(n+1)-time_vector(n)))/(latent_heat_evap_hydrogen*number_of_tanks*volume_tank_first)-storing_pressure
    pressure_gas_vap(n) = (R_constant*Temp_saturation(n)*Q_tot(n)*(time_vector(n+1)-time_vector(n)))/(latent_heat_evap_hydrogen*number_of_tanks*volume_tank_first);
    
    if pressure_gas_vap(n)>=2*storing_pressure_first
        storing_pressure = storing_pressure_first;
    else 
        storing_pressure = pressure_gas(n);
    end 
    volume_gas(n) = 
    
    BOM(n) = (pressure_gas(n)*number_of_tanks*volume_tank_first)/(R_constant*Temp_saturation(n));
end 

BOM_important = BOM(30:end);
time_start = time_vector(1);
new_time = [];
for j = 1:length(R)
    new_time(j) = time_start + time_vector(j);
end 
time_important = new_time(30:end);
Boil_off_mass = real(sum(BOM_important)); %kg
BOM_mt = Boil_off_mass/1000; %MT
end 

