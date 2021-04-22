%jjs280
%03/31/2021
%This code is a wrapper used to model and size the NEP system
close all
clear
clc

%Optimization Parameters are Produced Thrust (T) and Input Power (P)
%Variables of interest are voltage (V), current (I), and specific impulse (Isp)

%Electrical Power: P = Vb*Ib
%Thrust for Xenon Propellant (Goebel Equation 2.3-17): T = 1.65*gamma*Ib*sqrt(Vb);
%Isp for Xenon Propellant (Goebel Equation 2.4-10): Isp = 123.6*gamma*eta_m*sqrt(Vb)

%Generating throttle curves (Thrust vs. Power, Isp vs Power)
%Performance constants:
gamma = 0.958; %Thrust correction factor (10deg half angle beam, 10% single to double ionized ratio)
eta_m = 0.9; %Mass utilization factor
propellant = 'xenon';

%Define current and voltage inputs
Vb = 1.3e3:1e1:3e3; %[V]
Ib = 1:1:2e1; %[A]

%Generate throttle curves at various current across a voltage range
T = []; P = []; Isp = [];
for i = 1:1:length(Ib)
    T(:,i) = 1.65.*gamma.*Ib(i).*sqrt(Vb); %[mN]
    P(:,i) = Vb.*Ib(i); %[W]
    Isp(:,i) = 123.6.*gamma.*eta_m.*sqrt(Vb); %[s]
end

%Selected Operating Point (112 Thrusters)
V_op = 2250;
I_op = 12.3;
T_op = 1.65*gamma*I_op*sqrt(V_op); %[mN]
P_op = V_op*I_op; %[W]
Isp_op = 123.6*gamma*eta_m*sqrt(V_op); %[s]

disp(' ')
fprintf('Target Operating Point:\n')
fprintf('Required Power (W): %.3f\n',P_op)
fprintf('Produced Thrust (mN): %.3f\n',T_op)
num_thrusters = 112;
fprintf('Total Thrust (N): %.3f\n',T_op*num_thrusters/1000)
fprintf('Operating Isp (s): %.3f\n',Isp_op)

plot_throttle_curves = 0;

%% Plots
plot_results = 0;

%Plot results vs. Voltage and Current
if plot_results == 1
    figure()
    plot(Vb,P);
    grid on
    hold on
    title('Power at Fixed Current vs Voltage');
    xlabel('Volts');
    ylabel('Power [W]');
    %legend('');
    figure()
    plot(Vb,T);
    grid on
    hold on
    title('Thrust at Fixed Current vs Voltage');
    xlabel('Volts');
    ylabel('Thrust [mN]');
    %legend('');
    figure()
    plot(Vb,Isp);
    grid on
    hold on
    title('Isp at Fixed Current vs Voltage');
    xlabel('Volts');
    ylabel('Isp [s]');
    %legend('');
end

if plot_results == 1
    %Plot Throttle Curves
    figure()
    plot(P,T);
    grid on
    hold on
    title('Power and Thrust (Varying Voltage)');
    xlabel('Power [W]');
    ylabel('Thrust [mN]');
    xline(25000,'-','25kW');
    yline(1000,'-','1N');
    %legend('');
    figure()
    plot(P,Isp);
    grid on
    hold on
    title('Power and Isp (Varying Voltage)');
    xlabel('Power [W]');
    ylabel('Isp [s]');
    xline(25000,'-','25kW');
    yline(4000,'-','4000s');
    %legend('');
end
