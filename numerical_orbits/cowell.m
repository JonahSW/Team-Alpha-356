%jjs280
%03/09/2021
%Initial version of Cowell's method of orbital perturbations for low thrust trajectory calculations
%Uses no solver at the moment
clear
clc
close all

%Generate ephemeris data for planetary bodies
ephemerides = ephemerides(0);
%Solar System Constants
mu_sun = 2.9591220823e-4;%Solar Gravitational Parameter [AU^3/day^2]

%Define Initial Conditions
%Start near Earth orbit
r0 = ephemerides(6,1);
theta0 = ephemerides(5,1);
rdot0 = ephemerides(6,2)-ephemerides(6,1);
thetadot0 = ephemerides(5,2)-ephemerides(5,1);

%Define Simulation Parameters
%initial time [days]
t0 = 0;
%final time [days]
t1 = 366;
%timestep [days]
dt = 30*(1/86400);
%divergence condition
divtheta = 1e3;
divr = 0.1;

a = 0;%Perturbation is zero

%Numerical Integration of Cowell's Method
r_output = zeros((t1-t0), 2);%Output position vector
rdot_output = zeros((t1-t0), 2);%Output velocity vector
r = r0;%Initial condition
theta = theta0;%Initial condition
t = t0;%Initial condition
rdot = rdot0;%Initial condition
thetadot = thetadot0;%Initial condition

%Derivatives
rdoubledot = 0;
thetadoubledot = 0;

index = 0;
while t <= t1
    index = index+1;
    %Update Acceleration
    rdoubledot = (-mu_sun/r^2) + a;
    thetadoubledot = (mu_sun/r^3)*theta + a;
    %Check if acceleration is diverging
    if((thetadoubledot>divtheta) || (rdoubledot>divr))
        t = t1;
    else
        %Update Velocity
        rdot = rdot + dt*rdoubledot;
        thetadot = thetadot + dt*thetadoubledot;

        %Update Position
        r = r + dt*rdot;
        theta = theta + dt*thetadot;

        r_output(index,1) = theta;
        r_output(index,2) = r;
        rdot_output(index,1) = thetadot;
        rdot_output(index,2) = rdot;
        t = t+dt
    end
end

% Plot Results

figure(5)
polarplot(ephemerides(1,:),ephemerides(2,:),ephemerides(3,:),ephemerides(4,:),ephemerides(5,:),...
    ephemerides(6,:),ephemerides(7,:),ephemerides(8,:),ephemerides(9,:),ephemerides(10,:));
hold on
polarplot(r_output(:,1),r_output(:,2),'LineWidth',2)


