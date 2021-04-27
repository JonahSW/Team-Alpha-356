%jjs280
%04/13/2021
%This wrapper code loads a spacecraft input and then attempts to generate an optimized trajectory
%rk4 solver and basic equations come from low_thrust.m and rk4.m code written by Dr. Barnhart
%
close all
clear
clc

%% Start up and read from input file
disp('low_thrust_modified');
fnmin = input('enter input file -> ','s'); fin=fopen(fnmin);
write = input('Write final ephemeris to output file -> ','s');
if (strcmp(write,'yes')) || (strcmp(write,'Yes')) || (strcmp(write,'y')) || (strcmp(write,'Y'))
    fnmout=input('enter output file -> ','s'); fout=fopen(fnmout,'w');
    disp(' ');
end
fprintf('Inputs:\n')
%...read the input file, write to display
%Define simulation parameters
h=fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',h,s);
iprt=fscanf(fin,'%hd',[1,1]); s=fgetl(fin); fprintf(1,'%hd %s\n',iprt,s);
n=fscanf(fin,'%hd',[1,1]); s=fgetl(fin); fprintf(1,'%hd %s\n',n,s);
mu = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',mu,s);

%Define trajectory parameters
r0 = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',r0,s);
r_target = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',r_target,s);
delta_i = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',delta_i,s);
theta_dot_mod = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',theta_dot_mod,s);
r_dot_initial = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',r_dot_initial,s);
t_thrust = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',t_thrust,s);
t_coast = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',t_coast,s);

%Define spacecraft parameters
m = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',m,s);
thrust = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',thrust,s);
Isp = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',Isp,s);

%% Calculate Non-dimensional Parameters
%Time
x1 = 0;
x2 = t_thrust*24*3600*sqrt(mu/r0^3);
x3 = t_coast*24*3600*sqrt(mu/r0^3);
%Initial Thrust
nu = (thrust/m)*(r0^2/mu)*1e-3;% converting mu and r0 from km to m
%Radial and Tangential Velocities
Vcir_initial = sqrt(mu/r0);
theta_dot_initial = Vcir_initial/r0;
rho_prime_initial = (r_dot_initial/r0)*sqrt(r0^3/mu); % rho prime initial
theta_prime_initial = theta_dot_initial*sqrt(r0^3/mu) + (theta_dot_mod/r0)/sqrt(mu/r0^3); % theta prime initial including mod
rho_initial = 1;
theta_initial = 0;
y = zeros(n);
y(1) = rho_prime_initial;
y(2) = rho_initial;
y(3) = theta_prime_initial;
y(4) = theta_initial;
rho_target = r_target/r0;

% Display initial nondimensional parameters
display_initial_parameters = 1;
if display_initial_parameters == 1
    disp(' ')
    fprintf('Nondimensional parameters:\n')
    fprintf('Target non-dimensional radius (rho): %.3f \n',rho_target)
    fprintf('Non-Dimen. Thrust (nu): %.6f \n',nu)
    fprintf('Non-Dimen. Total Time (x2+x3): %.3f \n',(x2+x3))
    disp(' ')
end

%% 
