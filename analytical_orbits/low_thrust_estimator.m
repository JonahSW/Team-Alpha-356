%jjs280
%03/13/2021
%Code for estimating the properties of a low thrust trajectory
%Takes the deltaV, Isp, Thrust, and Mass Flow Rate of a low thrust trajectory and engine and
%estimates the trajectories duration

function duration = low_thrust_estimator(deltaV, Isp, Fthrust, MFR)
    duration = Isp*Fthrust/(deltaV*MFR);
end