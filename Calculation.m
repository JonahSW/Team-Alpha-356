function [BoiloffMasstoCeres] = Calculation(a_1,a_2,t_T,total_mass_fraction,M_d,r_tank,number_of_tanks)

% a_1 = 1; % Distance from Earth to Sun 
% a_2 = 2.766; % Distance from Ceres to Sun 
% t_T = 1.2919;
% total_mass_fraction = 2.7472;
% M_d = 167.056;
% r_tank = 4.4;
% number_of_tanks = 3;

a_u = 1.5E+8; % AU to Km
distance_C = a_2 * a_u; % in km
distance_C_m = distance_C*1000;
distance_E_m = a_1*a_u *1000;
R_s_p = distance_E_m:100000000:distance_C_m;
Temp_2 = 2.7; %k

[distance_from_sun, solar_flux] = Sun_radiation(a_1,a_2);

[M_w_first, h_tank_first, volume_tank_first, surface_area, thickness] = tank_size(total_mass_fraction, M_d, r_tank, number_of_tanks);

[R_s_p, R, q_tot] = heat_transfer(M_w_first, a_2, a_1, solar_flux, surface_area, thickness, h_tank_first, Temp_2);

[BOM_important, time_important, BOM_mt] = Boiloffmassguy(t_T, R_s_p, q_tot, volume_tank_first, R, number_of_tanks);

[R_S_P, BOM_Rate] = Boiloffmassrate(BOM_important, time_important, R_s_p);

BoiloffMasstoCeres = BOM_mt;


Distance = distance_from_sun((30:length(R_s_p)))/(a_u *1000);
figure(1)
yyaxis left
plot(time_important/3.154E7,BOM_important)
ylabel('BOM')
hold on 
yyaxis right
plot(time_important/3.154E7,Distance)
ylabel('Distance (m)')
xlabel('Time (years)')

figure(2)
plot(R_S_P,BOM_Rate)
xlabel('R_{sp}')
ylabel('BOM rate')
end 