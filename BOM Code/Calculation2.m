function [BoiloffMasstoEarth] = Calculation2(a_1,a_2,t_T,total_mass_fraction,M_d,r_tank,number_of_tanks)


% a_1 = 1; % Distance from Earth to Sun 
% a_2 = 2.766; % Distance from Ceres to Sun 
% t_T = 1.2919;
% total_mass_fraction = 2.7249;
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

[R_s_p, R, q_tot] = heat_transfer(M_w_first, a_2, a_1, flip(solar_flux), surface_area, thickness, h_tank_first, Temp_2);

[BOM_important, time_important, BOM_mt] = Boiloffmassguy(t_T, flip(R_s_p), q_tot, volume_tank_first, R, number_of_tanks);

[R_S_P, BOM_Rate] = Boiloffmassrate(BOM_important, time_important, flip(R_s_p));

BoiloffMasstoEarth = BOM_mt;

Distance = flip(distance_from_sun((30:length(R_s_p)))/(a_u *1000));
figure(1)
h1 = axes;
yyaxis left
plot(flip(time_important)/3.154E7,BOM_important)
ylabel('BOM')
hold on 
yyaxis right
plot(flip(time_important)/3.154E7,Distance)
ylabel('Distance (AU)')
xlabel('Time (years)')
set(h1, 'Xdir', 'reverse')

figure(2)
h2 = axes;
plot(flip(R_S_P),BOM_Rate)
xlabel('R_{sp}')
ylabel('BOM rate')
set(h2, 'Xdir', 'reverse')
