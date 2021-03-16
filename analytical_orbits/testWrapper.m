%jjs280
%03/14/2021
%Wrapper file for testing analytical orbit calculation tools
close all
clc
clear

%% Testing orbit calculations

%Load ephemeris data
ephemerides = ephemerides;

a_Earth = 1.496e+8*mean(ephemerides(4,:));%Convert AU to km
a_Ceres = 1.496e+8*mean(ephemerides(8,:));%Convert AU to km
    
Isp = 950;

mu_sun = 1.32712e11;%Solar Gravitational Parameter [km^3/sec^2]
L12_Hohmann = pi*(1-((a_Earth+a_Ceres)/(2*a_Ceres))^(3/2));%Determines starting separation angle for Hohmann transfer
thetaB_prime = pi;%Limiting case Hohmann Transfer

[Mratio_A, Mratio_B, deltaV_A, deltaV_B, tT] = elliptic_transfer(a_Earth, a_Ceres, thetaB_prime, Isp);


%% Testing low thrust estimator
deltaV = 9500;%[m/s]
Isp = 1000:1000:6000;%[s]
M_dry = 200000;%[kg]
thrust = 100e-3:10e-3:100;%100mN to 100N

duration1 = low_thrust_estimator(deltaV, Isp(1), M_dry, thrust);
duration2 = low_thrust_estimator(deltaV, Isp(2), M_dry, thrust);
duration3 = low_thrust_estimator(deltaV, Isp(3), M_dry, thrust);
duration4 = low_thrust_estimator(deltaV, Isp(4), M_dry, thrust);
duration5 = low_thrust_estimator(deltaV, Isp(5), M_dry, thrust);
duration6 = low_thrust_estimator(deltaV, Isp(6), M_dry, thrust);
duration = vertcat(duration1,duration2,duration3,duration4,duration5,duration6);

figure(1)
semilogy(thrust,duration);
grid minor
title({['Thrust vs. Trajectory Duration'],...
    ['(deltaV = ',num2str(deltaV),'m/s and Mdry = ',num2str(M_dry),' kg)']});
xlabel('Thrust (N)');
ylabel('Trajectory Duration (days)');
yline(365.25*5,'-','5 Years')
yline(365.25*2,'-','2 Years')
legend('Isp = 1000s','Isp = 2000s','Isp = 3000s','Isp =4000s','Isp = 5000s','Isp = 6000s');

