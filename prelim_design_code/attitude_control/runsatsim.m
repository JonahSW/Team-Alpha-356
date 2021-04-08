%% Run file for satsim.m
clear
clc
close all
t = 0:0.01:40;
Kp = 36;
Ki = 0;
Kd = 3.3;
gains = [Kp Ki Kd];

disp('Maximum Requirements:')
disp('    Overshoot Settling  Speed     Torque    Hp        BW        Energy');
disp('      0.20      20.0    500.0     6.565    0.625    ------      ------');

disp('Output:')
satsim(t,gains);