%jjs280
%03/14/2021
%Wrapper file for testing analytical orbit calculation tools
close all
clc
clear

%% Testing orbit calculations

%Load ephemeris data
positions = planetary_positions(1,0);
earth_ephemeris = read_ephemeris(2);
ceres_ephemeris = read_ephemeris(4);

%Calculate departure windows for Hohmann Transfers
[departure_days, departure_dates] = hohmann_window(2,4,1);

Isp = 950;%Isp for NTP

[Mratio_A, Mratio_B, deltaV_A, deltaV_B, tT] = elliptic_transfer(earth_ephemeris(3,912), ceres_ephemeris(3,912), pi, Isp)

%Calculate Mass Ratios and deltaVs for each launch window
half_window_width = 7;% Days on either side of launch window [days]
for i = 1:1:length(departure_days)
    
    for ii = 1:1:2*half_window_width
    %[Mratio_A, Mratio_B, deltaV_A, deltaV_B, tT] = elliptic_transfer(a_Earth, a_Ceres, thetaB_prime, Isp);
    %Mratio(i) = Mratio_A + Mratio_B;
    %deltaV(i) = deltaV_A + deltaV_B;
    end
    
end

%% Testing low thrust estimator
deltaV = 9750;%[m/s]
Isp = 1000:1000:6000;%[s]
M_dry = 100e3:20e3:200e3;%[kg]
thrust = 100e-3:10e-3:100;%100mN to 100N

%Varying Isp
duration1 = low_thrust_estimator(deltaV, Isp(1), M_dry(4), thrust);
duration2 = low_thrust_estimator(deltaV, Isp(2), M_dry(4), thrust);
duration3 = low_thrust_estimator(deltaV, Isp(3), M_dry(4), thrust);
duration4 = low_thrust_estimator(deltaV, Isp(4), M_dry(4), thrust);
duration5 = low_thrust_estimator(deltaV, Isp(5), M_dry(4), thrust);
duration6 = low_thrust_estimator(deltaV, Isp(6), M_dry(4), thrust);
duration = vertcat(duration1,duration2,duration3,duration4,duration5,duration6);

plot_lte = 1;
if plot_lte == 1
    figure()
    semilogy(thrust,duration);
    grid minor
    title({['Thrust vs. Trajectory Duration'],...
        ['(deltaV = ',num2str(deltaV),'m/s and Mdry = ',num2str(M_dry(4)),' kg)']});
    xlabel('Thrust (N)');
    ylabel('Trajectory Duration (days)');
    yline(365.25*5,'-','5 Years')
    yline(365.25*2,'-','2 Years')
    legend('Isp = 1000s','Isp = 2000s','Isp = 3000s','Isp =4000s','Isp = 5000s','Isp = 6000s');
end

%Varying dry mass
duration1 = low_thrust_estimator(deltaV, Isp(4), M_dry(1), thrust);
duration2 = low_thrust_estimator(deltaV, Isp(4), M_dry(2), thrust);
duration3 = low_thrust_estimator(deltaV, Isp(4), M_dry(3), thrust);
duration4 = low_thrust_estimator(deltaV, Isp(4), M_dry(4), thrust);
duration5 = low_thrust_estimator(deltaV, Isp(4), M_dry(5), thrust);
duration6 = low_thrust_estimator(deltaV, Isp(4), M_dry(6), thrust);
duration = vertcat(duration1,duration2,duration3,duration4,duration5,duration6);

if plot_lte == 1
    figure()
    semilogy(thrust,duration);
    grid minor
    title({['Dry Mass vs. Trajectory Duration'],...
        ['(deltaV = ',num2str(deltaV),'m/s and Isp = ',num2str(Isp(4)),' s)']});
    xlabel('Dry Mass (mT)');
    ylabel('Trajectory Duration (days)');
    yline(365.25*5,'-','5 Years')
    yline(365.25*2,'-','2 Years')
    legend('Mdry = 100mT','Mdry = 120mT','Mdry = 140mT','Mdry = 160mT','Mdry = 180mT','Mdry = 200mT');
end
