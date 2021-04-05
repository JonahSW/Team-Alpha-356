function [BOM_atCeres] = Calculation3(a_1,a_2,t_T,total_mass_fraction,M_d,r_tank,number_of_tanks,str)


% total_mass_fraction = 2.72
% M_d = 88.155 MT
% r_tank = 4.4
% number_of_tanks = 2
%a_1 = 2.766;
%a_2 = 2.766;
%t_T = 1.276;
a_u = 1.60E11; %AU to m
years_to_secs = 3.15E7;

%% Inputs from other files

[surface_area] = tank_size(total_mass_fraction, M_d, r_tank, number_of_tanks,str);
% total_mass_fraction = 2.72
% M_d = 88.155 MT
% r_tank = 4.4
% number_of_tanks = 2


[distance_from_sun,Boil_off_rate,BOM,BOM_duration] = BoiloffMasslockheed(str,a_1,a_2,t_T,surface_area);

disp('--------')
BOM_atCeres = BOM_duration;


time =  1:((t_T*years_to_secs)/length(distance_from_sun)):t_T*years_to_secs;
figure()
plot(distance_from_sun./a_u,Boil_off_rate)
title('Orbiting Ceres')
ylabel('Boil off Mass Rate (kg/s)')
xlabel('Distance (AU)')

figure()
plot(time./years_to_secs,BOM)
title('Orbiting Ceres')
ylabel('Boil off Mass (kg)')
xlabel('Time (Years)')
end
