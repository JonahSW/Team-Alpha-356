%--------------------------------------------------------------------------
% low_thrust_modified - low/continuous thrust orbit solution (using rk4)
% edited by sdf33, jjs280
%--------------------------------------------------------------------------
clear; close all; clc;
disp('low_thrust_modified');
fnmin=input('enter input file -> ','s'); fin=fopen(fnmin);
fnmout=strcat(fnmin(1:end-3),'.out'); %...input('enter output file -> ','s');
fout=fopen(fnmout,'w');
disp(' ');
fprintf('Inputs:\n')
%...read the input file, write to display
itan=fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',itan,s);
nu=fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',nu,s);
x1=fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',x1,s);
x2=fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',x2,s);
x3=fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',x3,s);
h=fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',h,s);
iprt=fscanf(fin,'%hd',[1,1]); s=fgetl(fin); fprintf(1,'%hd %s\n',iprt,s);
n=fscanf(fin,'%hd',[1,1]); s=fgetl(fin); fprintf(1,'%hd %s\n',n,s);

for i=1:n
    yv=fscanf(fin,'%g',[1,1]); s=fgetl(fin);
    fprintf(1,'%g %s\n',yv,s);
    y(i)=yv;
end

%%...begin edits....
r0 = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',r0,s);
m = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',m,s);
mu = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',mu,s);
r_target = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',r_target,s);
V_kick = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',V_kick,s);
r_dot_initial = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',r_dot_initial,s);
theta_dot_initial = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',theta_dot_initial,s);
t_thrust = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',t_thrust,s);
t_coast = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',t_coast,s);
thrust = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',thrust,s);
Isp = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',Isp,s);
plane_change = fscanf(fin,'%g',[1,1]); s=fgetl(fin); fprintf(1,'%g %s\n',plane_change,s);

%...Replace x2 and x3 with values calculated based on time for thrusting and coasting (Converting from days to seconds)
x2 = t_thrust*24*3600*sqrt(mu/r0^3);
x3 = t_coast*24*3600*sqrt(mu/r0^3);

%...Circular Orbit Velocity
Vcir = sqrt(mu/r0);
theta_dot_initial = Vcir/r0;

%...Replace nu with value calculated based maximum thrust
nu = (thrust/m)*(r0^2/mu)*1e-3;%...converting from km to m
rho_prime_initial = (r_dot_initial/r0)*sqrt(r0^3/mu); %...rho prime initial
theta_prime_initial = theta_dot_initial*sqrt(r0^3/mu); %...theta prime initial
y(1) = y(1) + rho_prime_initial; %...corrected rho prime
y(3) = y(3) + theta_prime_initial + ((V_kick)/r0)/sqrt(mu/r0^3); %...corrected theta prime
V0 = norm(r_dot_initial,r0*theta_dot_initial); %...initial velocity magnitude

%...Display some initial parameters
disp(' ')
disp(' ')
fprintf('Outputs:\n')
rho_target = r_target/r0;
fprintf('Target rho: %.3f \n',rho_target)
fprintf('Non-Dimen. Thrust (nu): %.6f \n',nu)
fprintf('Thrust: %.3f N\n',thrust)
fprintf('Non-Dimen. Total Time (x2+x3): %.3f \n',x2+x3)
fprintf('Total Time: %.3f days\n',t_thrust+t_coast)
fprintf('Initial Velocity w/ Vkick: %.3f km/s\n',V0+V_kick)

%% end edits......more edits below

x=x1; 
yrow=y(1:end); %fprintf(fout,'%12.4e',x); fprintf(fout,'%12.4e',yrow); 
%fprintf(fout,'\r\n');
ii=0;
nn=0;
rho(1)=y(2); theta(1)=y(4);
while x <= x2-h
    ii=ii+1;
    dydx=derivs(x,y,nu,itan);
    yout=rk4(y,dydx,n,x,h,nu,itan);
    if (ii==iprt)
        yrow=yout(1:end); xp=x+h; %fprintf(fout,'%12.4e',xp); 
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
%...write to display C3 and state variables at the end of thrusting
c3=y(1)*y(1)+y(2)*y(2)*y(3)*y(3)-2.0/y(2); yrow=y(1:end);
%fprintf(1,'\n');
%fprintf(1,'%12.4e',c3); fprintf(1,'%12.4e',xp); fprintf(1,'%12.4e',yrow);
%fprintf(1,'\r\n');
rho_end_thrust=y(2); theta_end_thrust=y(4);
%...zero-thrust coasting portion of the trajectory
nu=0.0;
while x < x2+x3
    ii=ii+1;
    dydx=derivs(x,y,nu,itan);
    yout=rk4(y,dydx,n,x,h,nu,itan);
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
rho_end=y(2); theta_end=y(4);
%...plot the spiral orbit in polar coordinates
%r1=linspace(1,1); th=linspace(0,2*pi);
%polarplot(th,r1,'-r',theta,rho,'-b',theta_end,rho_end,'ro');
%...if version of MATLAB does not support polarplot, use polar below
%polar(theta,rho);
status=fclose('all');


%%...begin edits......

%...Display some final results
r_final = r0*rho_end;
disp(' ')
fprintf('Initial Orbit Radius:   %.3f km\n',r0)
fprintf('Final Orbit Radius:     %.3f km\n',r_final)
fprintf('Desired Orbit Radius:   %.3f km\n',r_target)

%...final dimen derivatives
r_dot_final = r0*y(1)*sqrt(mu/r0^3); %...final dimensional radial velocity
r_final = r0*rho_end; %...final dimensional radius
theta_dot_final = y(3)*sqrt(mu/r0^3); %...final dimensional circumferential velocity
Vf = norm(r_dot_final,r_final*theta_dot_final); %...final velocity

fprintf('\nFinal Dimensional Derivatives:\n')
fprintf('r dot final:   %.5f km/s\n',r_dot_final)
fprintf('theta dot final:   %.5f rad/s\n',theta_dot_final)

%...final non-dimen derivatives
rho_prime_end = y(1); %...final non-dimen rho deriv
theta_end = y(4); %...final non-dimen theta
theta_prime_end = y(3); %...final non-dimen theta deriv

disp(' ')
fprintf('Final Non-Dimen Derivatives:\n')
fprintf('rho_end:         %.3f\n',rho_end)
fprintf('rho_prime_end:    %.3f\n',rho_prime_end)
fprintf('theta_end:       %.3f\n',theta_end)
fprintf('theta_prime_end:    %.3f\n',theta_prime_end)

C = y(1)^2 + (y(2)*y(3))^2 - 2/y(2); %...determine which final orbit type
disp(' ')
if C < 0
    fprintf('Final Orbit Type:   Elliptical\n')
elseif C == 0
    fprintf('Final Orbit Type:   Parabolic\n')
else
    fprintf('Final Orbit Type:   Hyperbolic\n')
end

fprintf('Final Velocity Magnitude:    %.3f km/s\n\n',Vf)


%...Calculate propellant mass consumed for burn:

g0 = 9.807; %[m/s^2]
propellant_mass = (thrust/(Isp*g0))*t_thrust*24*3600;
fprintf('Propellant Mass Consumed:   %.3f kg\n',abs(propellant_mass))
fprintf('Final Mass:   %.3f kg\n',m-abs(propellant_mass))

%%...plot the spiral orbit in polar coordinates
r1=linspace(1,1); th=linspace(0,2*pi); 
desired_r=linspace(r_target,r_target); desired_th=linspace(0,2*pi);
polarplot(th,r1,'-r',theta,rho,'-b',theta_end_thrust,rho_end_thrust,'ro',desired_th,desired_r/r0,'-g')
%%...end edits........
