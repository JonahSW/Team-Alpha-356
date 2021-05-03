%jjs280
%04/13/2021
%This wrapper code simulates the spacecraft's trajectory
close all
clear
clc

%% Constants

%Alpha Spacecraft
LEO_assembly_altitude = 800e3;% [m]
LEO_capture_altitude = 1000e3;% [m]
LCO_altitude = 700e3;%[m]
launch_inclination = 28.573*(pi/180);
Isp = 8426;% [s]
Isp_kick = 302.9;% [s]
g0 = 9.807;% [m/s^2]
m_dry = 300e3;% [kg]
m_dry_kickstage = 531.2;% [kg]
m_wet_kickstage = 5031.84;% [kg]
m_payload = 14636;% [kg]
m_sample = 50;% [kg]
m_dot_propellant = 1.51e-3;%[kg/s]
ceres_loiter_time = 60*(24*3600);% [s]
deltaV_kickstage = 100;% [m/s]
ep_thrust = 122.32;% [N]
num_boosters = 4;% [kg]

%Solar System
%Sun
mu_sun = 1.32712440018e20; %[m^3/s^2]
AU = 149597870.7e3;% [m]
%Earth
mu_earth = 3.986004418e14; %[m^3/s^2]
earth_ephemeris = read_ephemeris(2);
a_earth = mean(earth_ephemeris(3,:))*AU;
in_earth = mean(earth_ephemeris(4,:));
r_earth = 6378000;% [m]
r_earth_SOI = 0.929e9;% [m]
earth_axial_tilt = 23.43928*pi/180;% [rad]
%Ceres
mu_ceres = 6.26325e10; %[m^3/s^2]
r_ceres = 469730;% [m]
r_ceres_SOI = 7709455.15;% [m]
ceres_ephemeris = read_ephemeris(4);
a_ceres = mean(ceres_ephemeris(3,:))*AU;
in_ceres = mean(ceres_ephemeris(4,:));

%% Kick Stage Sizing (Partial kick)
%Chemical Kick Stage for LEO Departure
solid_propellant_mass = num_boosters*(m_wet_kickstage-m_dry_kickstage);% [kg]
LEO_departure_dry_mass = 620e3;% [kg] (UPDATE VALUES AS CALCS ITERATE)
deltaV_kick = Isp_kick*g0*log((LEO_departure_dry_mass+solid_propellant_mass)/LEO_departure_dry_mass);
disp(['Delta V for the kick stage is: ',num2str(deltaV_kick),' m/s']);

%% Estimate Delta Vs
%Assume a low thrust gradual spiral with a sub-optimal plane change

%DeltaV for leaving Earth's SOI (1)
r0 = LEO_assembly_altitude + r_earth;
r = r_earth_SOI;
v1 = sqrt(mu_earth/r0);
v2 = sqrt(mu_earth/r);
delta_i = (launch_inclination - earth_axial_tilt)*180/pi;% plane change for leaving Earth, deg
deltaV1 = spiral(mu_earth,r0,r) + inclination_change(v1,v2,delta_i);
%REPLACE WITH OUTBOUND CALC FROM LOW_THRUST_MODIFIED and ALPHA_EARTH_OUT.in
deltaV1 = 7122.479;% [m/s]
V1_excess = 809;% [m/s]

%DeltaV for traveling to Ceres (2)
r0 = a_earth;
r = a_ceres;
v1 = sqrt(mu_sun/r0);
v2 = sqrt(mu_sun/r);
delta_i = (in_earth + in_ceres)*180/pi;% plane around sun, deg
deltaV2 = spiral(mu_sun,r0,r) + inclination_change(v1,v2,delta_i) - V1_excess;

%DeltaV for capturing in Ceres SOI (3)
r0 = r_ceres_SOI;
r = LCO_altitude + r_ceres;
deltaV3 = spiral(mu_ceres,r0,r);

%DeltaV for leaving Ceres' SOI (4)
r0 = LCO_altitude + r_ceres;
r = r_ceres_SOI;
deltaV4 = spiral(mu_ceres,r0,r);
%REPLACE WITH OUTBOUND CALC FROM LOW_THRUST_MODIFIED and ALPHA_CERES_OUT.in
deltaV4 = 358.257;% [m/s]
V4_excess = 179;% [m/s]

%DeltaV for traveling to Earth (5)
r0 = a_ceres;
r = a_earth;
v1 = sqrt(mu_sun/r0);
v2 = sqrt(mu_sun/r);
delta_i = (in_earth + in_ceres)*180/pi;% plane around sun, deg
deltaV5 = spiral(mu_sun,r0,r) + inclination_change(v1,v2,delta_i) - V4_excess;

%DeltaV for capturing in Earth SOI (6)
r0 = r_earth_SOI;
r = LEO_capture_altitude + r_earth;
deltaV6 = spiral(mu_earth,r0,r);

%% Mission Phases
m_final = m_dry+m_sample-m_payload;% [kg]
%Earth Inbound
mratio_6 = exp(deltaV6/(Isp*g0));
m_wet_6 = m_final*mratio_6;% [kg]
m_propellant_6 = m_wet_6 - m_final;% [kg]
duration_6 = m_propellant_6/m_dot_propellant;% [s]
%Sun Inbound
mratio_5 = exp(deltaV5/(Isp*g0));
m_wet_5 = (m_wet_6)*mratio_5;% [kg]
m_propellant_5 = m_wet_5 - m_wet_6;% [kg]
duration_5 = m_propellant_5/m_dot_propellant;% [s]
%Ceres Outbound
mratio_4 = exp(deltaV4/(Isp*g0));
m_wet_4 = (m_wet_5)*mratio_4;% [kg]
m_propellant_4 = m_wet_4 - m_wet_5;% [kg]
duration_4 = m_propellant_4/m_dot_propellant;% [s]
%Ceres Orbit
m_wet_ceres = m_wet_4 + m_payload - m_sample;% [kg]
duration_ceres = ceres_loiter_time;% [s]
%Ceres Inbound
mratio_3 = exp(deltaV3/(Isp*g0));
m_wet_3 = (m_wet_ceres)*mratio_3;% [kg]
m_propellant_3 = m_wet_3 - m_wet_ceres;% [kg]
duration_3 = m_propellant_3/m_dot_propellant;% [s]
%Sun Outbound
mratio_2 = exp(deltaV2/(Isp*g0));
m_wet_2 = (m_wet_3)*mratio_2;% [kg]
m_propellant_2 = m_wet_2 - m_wet_3;% [kg]
duration_2 = m_propellant_2/m_dot_propellant;% [s]
%Earth Outbound
mratio_1 = exp(deltaV1/(Isp*g0));
m_wet_1 = (m_wet_2)*mratio_1;% [kg]
m_propellant_1 = m_wet_1 - m_wet_1;% [kg]
duration_1 = m_propellant_1/m_dot_propellant;% [s]

m_wet_initial = m_wet_1 + m_wet_kickstage;% [kg]
%mratio1 = exp(deltaV1/(Isp*g0))
%% Results
deltaV_total = deltaV1+deltaV2+deltaV3+deltaV4+deltaV5+deltaV6;% [m/s]
mratio_total = mratio_1*mratio_2*mratio_3*mratio_4*mratio_5*mratio_6;
m_propellant_total = (m_propellant_1+m_propellant_2+m_propellant_3+m_propellant_4+m_propellant_5+m_propellant_6);% [kg]
duration_thrust = (duration_1+duration_2+duration_3+duration_4+duration_5+duration_6)/(24*3600);% [days] 
duration_total = (duration_1+duration_2+duration_3+duration_ceres+duration_4+duration_5+duration_6)/(24*3600);% [days]
disp(['The total mission deltaV is ',num2str(deltaV_total*1e-3),' km/s.']);
disp(['The total mission mass ratio is ',num2str(mratio_total),'']);
disp(['The total spacecraft wet mass is ',num2str(m_wet_1),' kg.']);
disp(['The total propellant mass is ',num2str(m_propellant_total),' kg.']);
disp(['The total thrusting duration is ',num2str(duration_thrust),' days.']);
disp(['The total mission duration is ',num2str(duration_total),' days.']);

figure()
durations = [0,duration_1,duration_2,duration_3,duration_ceres,duration_4,duration_5,duration_6];
masses = [m_wet_initial,m_wet_1,m_wet_2,m_wet_3,m_wet_ceres,m_wet_4,m_wet_5,m_wet_6,m_dry];
plot(1:1:length(masses),masses,'o');
grid minor
xline(1,'-','LEO Assembly Mass');
xline(2,'-','After Kickstage Burn','LabelVerticalAlignment','top');
xline(3,'-','After LEO Spiral Out');
xline(4,'-','After Sun Spiral Out');
xline(5,'-','After Ceres Spiral In');
xline(6,'-','At Ceres Departure');
xline(7,'-','After Ceres Spiral Out');
xline(8,'-','After Sun Spiral In');
xline(9,'-','After LEO Spiral In','LabelHorizontalAlignment','left');
title('Mass over Mission Duration');
ylabel('Mass (kg)');
xlabel('Time');

%% Functions

%Calculates the deltaV for a spiral orbit from r0 to r about a body with standard gravitational parameter mu
function deltaV = spiral(mu,r0,r)
    deltaV = abs(sqrt(mu/r0) - sqrt(mu/r));
end

%Calculates the deltaV for a sub optimal low thrust plane change given the initial and final velocities
function deltaV = inclination_change(v1,v2,delta_i)
    deltaV = sqrt(v1^2 + v2^2 - 2*v1*v2*cos((delta_i*pi/180)*pi/2));
end

%Implements a function that calculates the wait time for a given difference in true anomaly
%a1, a2, are semi major axes of initial and final heliocentric orbits
function [t_wait] = wait_time(a1, a2, dtheta)
    mu_sun = 1.32712e20;% Solar Gravitational Parameter [m^3/sec^2]
    %Convert to m from AU
    AU = 1.496e11;%AU in m
    a1 = a1*AU;
    a2 = a2*AU;
    %Calculations
    alpha = dtheta*(((a1+a2)/(2*a1))^(3/2)-1);
    S = (((2*pi)*sqrt(a1^3/mu_sun))^(-1) - ((2*pi)*sqrt(a2^3/mu_sun))^(-1))^(-1);%Synodic Period
    t_wait = (2-(alpha/pi))*(S/(24*3600));
end

%Implements a function that can calculate the departure date window an orbital transfer given the
%required difference in longitudes
%Returns an array of possible departure days as both indexes (days from 01/01/2040) and dates
%1 is Venus, 2 is Earth, 3 is Mars, 4 is Ceres, 5 in Jupiter
function [departure_days, departure_dates] = transfer_window(body1,body2,L12,plotme)

    body1_data = read_ephemeris(body1);
    body2_data = read_ephemeris(body2);
    
    a1 = body1_data(3,:);
    a2 = body2_data(3,:);
    theta1 = body1_data(2,:);
    theta2 = body2_data(2,:);
    
    departure_days = zeros(1,1);
    dtheta = zeros(1,1);
    num_windows = 1;
    for i = 1:1:length(a1)
        dtheta(i) = theta2(i) - theta1(i);
        if (dtheta(i) > (L12-L12*0.01)) && (dtheta(i) < (L12+L12*0.01))
            departure_days(num_windows) = i;
            %disp('Found One!');
            num_windows = num_windows+1;
        end
    end
    
    %Find time as dates
    t1 = datetime(2040,1,1,0,0,0);
    for i = 1:1:length(departure_days)
        departure_dates(i) = t1 + days(departure_days(i));
    end
    
    if plotme == 1
        figure()
        plot(1:1:length(theta1),theta1, 1:1:length(theta2),theta2);
        title('True Anomaly Over Time');
        xlabel('Time (days)');
        ylabel('True Anomaly (rad)');
        legend('Body 1','Body 2')
        
        figure()
        plot(1:1:length(dtheta),dtheta);
        hold on
        plot(1:1:length(L12),L12);
        title('Difference in True Anomaly Over Time');
        xlabel('Time (days)');
        legend('Difference in True Anomaly','Hohmann Transfer Alignment')
        
        for i = 1:1:length(departure_dates)
           xline(departure_days(i),'-',datestr(departure_dates(i)),'HandleVisibility','off'); 
        end
    end
end