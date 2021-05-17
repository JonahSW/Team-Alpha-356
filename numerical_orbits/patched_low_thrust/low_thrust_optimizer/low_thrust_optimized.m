%jjs280
%04/13/2021
%This wrapper code loads a spacecraft input and then attempts to generate an optimized trajectory
%rk4 solver and basic equations come from low_thrust.m and rk4.m code written by Dr. Barnhart
%
close all
clear
clc

%% Start up and read from input file
disp('low_thrust_optimized');
fnmin = input('enter input file -> ','s'); fin=fopen(fnmin);
write = input('Write final ephemeris to output file -> ','s');
if (strcmp(write,'yes')) || (strcmp(write,'Yes')) || (strcmp(write,'y')) || (strcmp(write,'Y'))
    fnmout=input('enter output file -> ','s'); fout=fopen(fnmout,'w');
    write = true;
    disp(' ');
else
    write = false;
end
fprintf('Inputs:\n')
%...read the input file, write to display
%Define simulation parameters
h=fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',h,s);
iprt=fscanf(fin,'%hd',[1,1]); s=fgetl(fin); fprintf(1,'%hd %s\n',iprt,s);
n=fscanf(fin,'%hd',[1,1]); s=fgetl(fin); fprintf(1,'%hd %s\n',n,s);
mu = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',mu,s);
final_orbit_type = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',final_orbit_type,s);
itan_initial = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',itan_initial,s);
max_iterations = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',max_iterations,s);

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

%% Set Optimizer
%End Condition: rho_end_coast = rho_target, rho_prime_final = 0;
%Define ending conditions
if final_orbit_type == 1
    rho_prime_final_target = 0;% For circular orbit
end

rho_prime_error = 0.05;
rho_prime_final = 1;
a = 0;
b = 0;
c = 1;
da = 1e-3;
db = 1e-3;
x = 0;

%Optimizer Loop
while (rho_prime_final < rho_prime_final_target*(1-rho_prime_error)) || (rho_prime_final > rho_prime_final_target*(1+rho_prime_error))
a = a+da;
b = b+db;
c = 0;
    
    %% Propagate Thrusting Trajectory
    rho_current = 0;
    x=x1; 
    yrow=y(1:end);
    if write
        fprintf(fout,'%12.4e',x); fprintf(fout,'%12.4e',yrow); fprintf(fout,'\r\n');
    end
    ii=0;
    nn=0;
    rho(1)=y(2); theta(1)=y(4);
    itan = itan_initial;
    while rho_current <= rho_target
        ii=ii+1;
        %Update itan according to steering law
        itan = a*y(2)^2 + b*y(2) + c;
        dydx=derivs(x,y,nu,itan);
        yout=rk4(y,dydx,n,x,h,nu,itan);
        if (ii==iprt)
            yrow=yout(1:end); xp=x+h;
            if write
                fprintf(fout,'%12.4e',xp); fprintf(fout,'%12.4e',yrow); fprintf(fout,'\r\n');
            end
            ii=0;
        end
    %...march the solution forward by interval h
        nn=nn+1;
        x=x1+h*nn;
        for i=1:n
            y(i)=yout(i);
        end
        rho(nn+1)=y(2); theta(nn+1)=y(4);
        rho_current=y(2);
    end
    %...write to display C3 and state variables at the end of thrusting
    c3=y(1)*y(1)+y(2)*y(2)*y(3)*y(3)-2.0/y(2); yrow=y(1:end);

    % Plot results of current trajectory
    figure(1);
    r1=linspace(1,1); th=linspace(0,2*pi); 
    desired_r=linspace(r_target,r_target); desired_th=linspace(0,2*pi);
    polarplot(th,r1,'-r',theta,rho,'-b',desired_th,desired_r/r0,'-g')
    
disp(a);disp(b);disp(c);    
% End Optimizer Loop
end

fprintf(1,'\n');
%fprintf(1,'%12.4e',c3); fprintf(1,'%12.4e',xp); fprintf(1,'%12.4e',yrow); fprintf(1,'\r\n');
rho_prime_final = y(1); rho_end_thrust=y(2); theta_end_thrust=y(4);


%% Propagate Final Orbit
 %...zero-thrust coasting portion of the trajectory
    nu=0.0;
    ii=0;
    nn=0;
    while x < x2+x3
        ii=ii+1;
        dydx=derivs(x,y,nu,0);
        yout=rk4(y,dydx,n,x,h,nu,0);
        if (ii==iprt)
            yrow=y(1:end); xp=x+h; %fprintf(fout,'%12.4e',xp); 
            %fprintf(fout,'%12.4e',yrow); fprintf(fout,'\r\n');
            ii=0;
        end
    %...march the solution forward by interval h
        nn=nn+1;
        x=x1+h*nn;
        for i=1:n
            y(i)=yout(i);
        end
        rho(nn+1)=y(2); theta(nn+1)=y(4);
    end
    rho_end_coast=y(2); theta_end_coast=y(4);

%% Plot Final Results
%%...plot the spiral orbit in polar coordinates
figure();
r1=linspace(1,1); th=linspace(0,2*pi); 
desired_r=linspace(r_target,r_target); desired_th=linspace(0,2*pi);
polarplot(th,r1,'-r',theta,rho,'-b',theta_end_thrust,rho_end_thrust,'ro',desired_th,desired_r/r0,'-g',theta_end_coast,rho_end_coast,'co')

%Print Key Results

%% Functions