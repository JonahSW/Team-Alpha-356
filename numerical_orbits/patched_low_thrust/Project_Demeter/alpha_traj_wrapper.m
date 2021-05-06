%jjs280
%04/13/2021
%This wrapper code simulates the spacecraft's trajectory
close all
clear
clc

%% Constants

%Alpha Spacecraft
Isp = 8372.1;% [s] (UPDATE REGULARLY)
Isp_kick = 302.9;% [s]
m_dot_propellant = 1.51e-3;%[kg/s]
ep_thrust = 122.32;% [N]
m_dry = 282e3;% [kg] (UPDATE REGULARLY) (Not including kickstages)
m_dry_kickstage = 507.16;% [kg]
m_wet_kickstage = 6742.16;% [kg]
m_payload = 15197.6;% [kg]
m_sample = 50;% [kg]
num_boosters = 4;% [kg]

LEO_departure_dry_mass = 437959.5897;% [kg] (UPDATE VALUES AS CALCS ITERATE)

LEO_assembly_altitude = 1100e3;% [m]
LEO_capture_altitude = 1000e3;% [m] UNUSED FOR LUNAR RENDEZVOUS
LLO_capture_altitude = 10000e3;% [m]
LCO_altitude = 700e3;%[m]
launch_inclination = 28.573*(pi/180);
ceres_loiter_time = 60;% [days]

LEO_assembly_duration = 30*6;% [days]
crew_earth_return_duration = 14;% [days]

%Solar System
%Sun
mu_sun = 1.32712440018e20; %[m^3/s^2]
AU = 149597870.7e3;% [m]
%Earth
g0 = 9.807;% [m/s^2]
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
%Moon
mu_moon = 4904869500000;% [m^3/s^2]
r_moon = 1737400;% [m]
r_moon_SOI = 66100000;% [m]
a_moon = 384399000;% [m]
in_moon = 5.145;% [deg]

%% Kick Stage Sizing (Partial kick)
%Chemical Kick Stage for LEO Departure
solid_propellant_mass = num_boosters*(m_wet_kickstage-m_dry_kickstage);% [kg]
deltaV_kick = Isp_kick*g0*log((LEO_departure_dry_mass+solid_propellant_mass)/LEO_departure_dry_mass);

%% Estimate Delta Vs
%Assume a low thrust gradual spiral with a sub-optimal climb and plane change

%DeltaV for leaving Earth's SOI (1)
r0 = LEO_assembly_altitude + r_earth;
r = r_earth_SOI;
delta_i = (launch_inclination - earth_axial_tilt)*180/pi;% plane change for leaving Earth, deg
deltaV1 = spiral_with_inclination_change(mu_earth,r0,r,delta_i) - deltaV_kick;
%REPLACE WITH OUTBOUND CALC FROM LOW_THRUST_MODIFIED and ALPHA_EARTH_OUT.in
%deltaV1 = 7122.479;% [m/s]
%V1_excess = 809;% [m/s]

%DeltaV for traveling to Ceres (2)
r0 = a_earth;
r = a_ceres;
delta_i = (in_earth + in_ceres)*180/pi;% plane change around sun, deg
deltaV2 = spiral_with_inclination_change(mu_sun,r0,r,delta_i);% - V1_excess;

%DeltaV for capturing in Ceres SOI (3)
r0 = r_ceres_SOI;
r = LCO_altitude + r_ceres;
deltaV3 = spiral(mu_ceres,r0,r);

%DeltaV for leaving Ceres' SOI (4)
r0 = LCO_altitude + r_ceres;
r = r_ceres_SOI;
deltaV4 = spiral(mu_ceres,r0,r);
%REPLACE WITH OUTBOUND CALC FROM LOW_THRUST_MODIFIED and ALPHA_CERES_OUT.in
%deltaV4 = 358.257;% [m/s]
%V4_excess = 179;% [m/s]

%DeltaV for traveling to Earth (5)
r0 = a_ceres;
r = a_earth;
delta_i = (in_earth + in_ceres)*180/pi;% plane change around sun, deg
deltaV5 = spiral_with_inclination_change(mu_sun,r0,r,delta_i);% - V4_excess;

%DeltaV for capturing in Earth SOI (6)
r0 = r_earth_SOI;
r = a_moon;
deltaV6 = spiral(mu_earth,r0,r);

%DeltaV for capturing in Lunar SOI (7)
r0 = a_moon;
r = LLO_capture_altitude;
delta_i = in_moon*pi/180;%Plane change to lunar orbit inclination, deg
deltaV7 = spiral_with_inclination_change(mu_moon,r0,r,delta_i);

%% Mission Phases
m_final = m_dry+m_sample-m_payload;% [kg]
%Moon Inbound
mratio_7 = exp(deltaV7/(Isp*g0));
m_wet_7 = m_final*mratio_7;% [kg]
m_propellant_7 = m_wet_7 - m_final;% [kg]
duration_7 = m_propellant_7/m_dot_propellant;% [s]
%Earth Inbound
mratio_6 = exp(deltaV6/(Isp*g0));
m_wet_6 = m_wet_7*mratio_6;% [kg]
m_propellant_6 = m_wet_6 - m_wet_7;% [kg]
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
m_propellant_1 = m_wet_1 - m_wet_2;% [kg]
duration_1 = m_propellant_1/m_dot_propellant;% [s]

%Assembly
m_wet_initial = m_wet_1 + m_wet_kickstage;% [kg]

%Overwride duration_2 (Sun outbound) and duration_5 (Sun inbound) from Zola Lengths
duration_2_acc = 362*24*3600;% [s]
duration_2_coast = 50*24*3600;% [s]
duration_2_dec = 80*24*3600;% [s]
duration_2 = duration_2_acc+duration_2_coast+duration_2_dec;

duration_5_acc = 235*24*3600;% [s]
duration_5_coast = 50*24*3600;% [s]
duration_5_dec = 80*24*3600;% [s]
duration_5 = duration_2_acc+duration_2_coast+duration_2_dec;

%Calculate Departure Date
dtheta_out_min = 281.064*pi/180;%Difference in true anomalies [rad]
dtheta_in_min = 97.1*pi/180;%Difference in true anomalies [rad]
dtheta_out_max = 310*pi/180;%Difference in true anomalies [rad]
dtheta_in_max = 107*pi/180;%Difference in true anomalies [rad]
t_wait = wait_time(a_earth, a_ceres, dtheta_in_min);
departure_dates = transfer_window(2,4,dtheta_out_min,dtheta_out_max,duration_2,1);

%% Results
deltaV_total = deltaV1+deltaV2+deltaV3+deltaV4+deltaV5+deltaV6+deltaV7;% [m/s]
mratio_total = mratio_1*mratio_2*mratio_3*mratio_4*mratio_5*mratio_6*mratio_7;
mratio_outbound = mratio_1*mratio_2*mratio_3;
mratio_inbound = mratio_4*mratio_5*mratio_6*mratio_7;
m_propellant_total = (m_propellant_1+m_propellant_2+m_propellant_3+m_propellant_4+m_propellant_5+m_propellant_6+m_propellant_7);% [kg]
duration_thrust = (duration_1+duration_2_acc+duration_2_dec+duration_3+duration_4+duration_5_acc+duration_5_dec+duration_6+duration_7)/(24*3600);% [days] 
duration_total = (LEO_assembly_duration+duration_thrust+duration_2_coast/(24*3600)+duration_5_coast/(24*3600)+duration_ceres+crew_earth_return_duration);% [days]
disp('Delta Vs:');
disp(['Delta V for the kick stage is: ',num2str(deltaV_kick),' m/s']);
disp(['The total mission deltaV is ',num2str(deltaV_total*1e-3),' km/s.']);
disp(['The outbound trajectory deltaV is ',num2str((deltaV1+deltaV2+deltaV3)*1e-3),' km/s.']);
disp(['The inbound trajectory deltaV is ',num2str((deltaV4+deltaV5+deltaV6+deltaV7)*1e-3),' km/s.']);
disp('Mass Ratios:');
disp(['The total mission mass ratio is ',num2str(mratio_total),'']);
disp(['The outbound trajectory mass ratio is ',num2str(mratio_outbound),'']);
disp(['The inbound trajectory mass ratio is ',num2str(mratio_inbound),'']);
disp('Masses:');
disp(['The total spacecraft wet mass is ',num2str(m_wet_1),' kg.']);
disp(['The total propellant mass is ',num2str(m_propellant_total),' kg.']);
disp('Durations:');
disp(['The total thrusting duration is ',num2str(duration_thrust),' days.']);
disp(['The total crewed mission duration is ',num2str(duration_total - LEO_assembly_duration),' days.']);
disp(['The total mission duration is ',num2str(duration_total),' days.']);

figure()
wet_masses = [m_wet_initial,m_wet_1,m_wet_2,m_wet_3,m_wet_ceres,m_wet_4,m_wet_5,m_wet_6,m_wet_7,m_final];
dry_masses = [m_dry+m_dry_kickstage,m_dry,m_dry,m_dry,m_dry-m_payload,m_dry-m_payload+m_sample,m_dry-m_payload+m_sample,m_dry-m_payload+m_sample,m_dry-m_payload+m_sample,m_dry-m_payload+m_sample];
plot(1:1:length(wet_masses),wet_masses,'-o',1:1:length(dry_masses),dry_masses,'-o');
grid minor
xline(1,'-','LEO Assembly Mass');
xline(2,'-','After Kickstage Burn','LabelVerticalAlignment','top');
xline(3,'-','After LEO Spiral Out');
xline(4,'-','After Sun Spiral Out');
xline(5,'-','After Ceres Spiral In');
xline(6,'-','At Ceres Departure');
xline(7,'-','After Ceres Spiral Out');
xline(8,'-','After Sun Spiral In');
xline(9,'-','After Earth Spiral In');
xline(10,'-','After Moon Spiral In','LabelHorizontalAlignment','left');
title('Mass over Mission Duration');
ylabel('Mass (kg)');
xlabel('Time');

figure()
durations = [LEO_assembly_duration,duration_1/(24*3600),duration_2/(24*3600),duration_3/(24*3600),...
    duration_ceres,duration_4/(24*3600),duration_5/(24*3600),duration_6/(24*3600),duration_7/(24*3600),crew_earth_return_duration];
plot(1:1:length(durations),durations,'o');
grid minor
xline(1,'-','LEO Assembly');
xline(3,'-','LEO Spiral Out');
xline(4,'-','Sun Spiral Out');
xline(5,'-','Ceres Spiral In');
xline(6,'-','Ceres Departure');
xline(7,'-','Ceres Spiral Out');
xline(8,'-','Spiral In');
xline(9,'-','Earth Spiral In');
xline(10,'-','Moon Spiral In','LabelHorizontalAlignment','left');
title('Mission Phase Duration');
ylabel('Duration (Days)');
xlabel('Time');

figure()
deltaVs = [deltaV1,deltaV2,deltaV3,deltaV4,deltaV5,deltaV6,deltaV7];
plot(1:1:length(deltaVs),deltaVs,'*r','linewidth',5);
grid minor
xline(1,'-','LEO Spiral Out');
xline(2,'-','Heliocentric Fast Spiral Out');
xline(3,'-','Ceres Spiral In');
xline(4,'-','Ceres Spiral Out');
xline(5,'-','Heliocentric Fast Spiral In');
xline(6,'-','Earth Spiral In');
xline(7,'-','Moon Spiral In','LabelHorizontalAlignment','left');
title('Delta V for Trajectory Segments');
ylabel('Delta V (m/s)');
xlabel('Trajectory Segment');
%% Functions

%Calculates the deltaV for a spiral orbit from r0 to r about a body with standard gravitational parameter mu
function deltaV = spiral(mu,r0,r)
    deltaV = abs(sqrt(mu/r0) - sqrt(mu/r));
end

%Calculates the deltaV for a sub optimal low thrust plane change given the initial and final velocities
function deltaV = spiral_with_inclination_change(mu,r0,r,delta_i)
    v1 = sqrt(mu/r0);
    v2 = sqrt(mu/r);
    deltaV = sqrt(v1^2 + v2^2 - 2*v1*v2*cos((delta_i*pi/180)*pi/2));
end

%Implements a function that calculates the wait time for a given difference in true anomaly
%a1, a2, are semi major axes of initial and final heliocentric orbits
function [t_wait] = wait_time(a1, a2, dtheta)
    mu_sun = 1.32712e20;% Solar Gravitational Parameter [m^3/sec^2]\
    %Calculations
    alpha = dtheta*(((a1+a2)/(2*a1))^(3/2)-1);
    S = (((2*pi)*sqrt(a1^3/mu_sun))^(-1) - ((2*pi)*sqrt(a2^3/mu_sun))^(-1))^(-1);%Synodic Period
    t_wait = (2-(alpha/pi))*(S/(24*3600));% [days]
end

%Implements a function that can calculate the departure date window an orbital transfer given the
%required difference in longitudes
%Returns an array of possible departure days as both indexes (days from 01/01/2040) and dates
%1 is Venus, 2 is Earth, 3 is Mars, 4 is Ceres, 5 in Jupiter
function [departure_days, departure_dates] = transfer_window(body1,body2,L12_min,L12_max,t_travel,plotme)

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
        if (dtheta(i) > (L12_min)) && (dtheta(i) < (L12_max))
            departure_days(num_windows) = i-t_travel;
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
        yline(L12_min,'-','Low Thrust Transfer Alignment','LabelVerticalAlignment','bottom')
        yline(L12_max,'-','Low Thrust Transfer Alignment')
        title('Difference in True Anomaly Over Time');
        xlabel('Time (days)');
        %legend('Difference in True Anomaly','Hohmann Transfer Alignment')
        
        xline(departure_days(1),'-',datestr(departure_dates(1)),'HandleVisibility','off','LabelHorizontalAlignment','left'); 
        xline(departure_days(length(departure_days)),'-',datestr(departure_dates(length(departure_days))),'HandleVisibility','off'); 
        
    end
end