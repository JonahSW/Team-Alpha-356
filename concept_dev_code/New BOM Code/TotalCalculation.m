function [total_BOM] = TotalCalculation()

str = input('What Element analysed for Boil Off: ','s')

%[BOM_toCeres] = Calculation(1,2.76,1.294,2.77,88.155,4.4,2,str)

%Hydrogen
% M_d = 88.155 MT
% r_tank = 4.4
% number_of_tanks = 2
%a_1 = 1;
%a_2 = 2.766;
%t_T = 1.294;
% total_mass_fraction = 2.77
% M_d = 88.155 MT
% r_tank = 4.4
% number_of_tanks = 2

%=====================
[BOM_toCeres] = Calculation(1,2.76,2.191,1.88,205.357,2,3,str);
%xenon
% M_d = 205.357 MT
% r_tank = 5
% number_of_tanks = 3
%a_1 = 1;
%a_2 = 2.766;
%t_T = 2.191;
% total_mass_fraction = 1.28
% number_of_tanks = 3

%[BOM_atEarth] = Calculation2(1,1,1/12,2.77,88.155,4.4,2,str)
% M_d = 88.155 MT
% r_tank = 4.4
% number_of_tanks = 2
%a_1 = 1;
%a_2 = 1;
%t_T = 1/12;
% total_mass_fraction = 2.77
% M_d = 88.155 MT
% r_tank = 4.4
% number_of_tanks = 2

%=====================
[BOM_atEarth] = Calculation2(1,1,2/12,1.88,205.357,2,3,str);

%xenon
% M_d = 205.357 MT
% r_tank = 5
% number_of_tanks = 3
%a_1 = 1;
%a_2 = 1;
%t_T = 2/12;
% total_mass_fraction = 1.28
% number_of_tanks = 3

%[BOM_atCeres] = Calculation3(2.766,2.766,.525,2.72,88.155,4.4,2,str)
% total_mass_fraction = 2.72
% M_d = 88.155 MT
% r_tank = 4.4
% number_of_tanks = 2
%a_1 = 2.766;
%a_2 = 2.766;
%t_T = .525;
% M_d = 88.155 MT
% r_tank = 4.4
% number_of_tanks = 2

%=====================
[BOM_atCeres] = Calculation3(2.766,2.766,0.493,1.88,205.357,2,3,str);
%xenon
% M_d = 205.357 MT
% r_tank = 5
% number_of_tanks = 3
%a_1 = 2.766;
%a_2 = 2.766;
%t_T = .493;
% total_mass_fraction = 1.28
% number_of_tanks = 3


%[BOM_toEarth] = Calculation4(1,2.766,1.294,2.72,88.155,4.4,2,str)
% total_mass_fraction = 2.72
% M_d = 88.155 MT
% r_tank = 4.4
% number_of_tanks = 2
%a_1 = 1;
%a_2 = 2.766;
%t_T = 1.294;
% M_d = 88.155 MT
% r_tank = 4.4
% number_of_tanks = 2

%======
[BOM_toEarth] = Calculation4(1,2.766,1.701,1.88,205.357,2,3,str);

%xenon
% M_d = 205.357 MT
% r_tank = 5
% number_of_tanks = 3
%a_1 = 1;
%a_2 = 2.766;
%t_T = 1.701;
% total_mass_fraction = 1.28
% number_of_tanks = 3

total_BOM = BOM_toCeres + BOM_atEarth + BOM_atCeres + BOM_toEarth;

end 