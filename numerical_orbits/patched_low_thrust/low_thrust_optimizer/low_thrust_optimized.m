%jjs280
%04/13/2021
%This wrapper code loads a spacecraft input and then attempts to generate an optimized trajectory
%rk4 solver and basic equations come from low_thrust.m and rk4.m code written by Dr. Barnhart
close all
clear
clc

%% Start up and read from input file
disp('low_thrust_modified');
fnmin=input('enter input file -> ','s'); fin=fopen(fnmin);
fnmout=strcat(fnmin(1:end-3),'.out'); %...input('enter output file -> ','s');
fout=fopen(fnmout,'w');
disp(' ');
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
theta_dot_initial = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',theta_dot_initial,s);
r_dot_initial = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',r_dot_initial,s);
t_thrust = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',t_thrust,s);
t_coast = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',t_coast,s);

delta_i = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',delta_i,s);
%Define spacecraft parameters
thrust = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',thrust,s);
Isp = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',Isp,s);
m = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',m,s);