%jjs280
%03/13/2021
%Code for estimating the properties of a low thrust trajectory
%Takes the deltaV [m/s], Isp [s], dry mass [kg], and Thrust [N] of a low thrust trajectory and engine and
%estimates the trajectories duration in days

function duration = low_thrust_estimator(deltaV, Isp, M_dry, thrust)
    g0 = 9.087;% [m/s^2]
    MFR = thrust./(g0*Isp);
    g0 = 9.807; %Std. gravitational acceleration [m/s^2]
    M_ratio = exp(deltaV./(Isp.*g0));
    M_wet = M_ratio.*M_dry;
    M_propellant = M_wet - M_dry;
    
    duration = (M_propellant./MFR)./(24.*3600);
end