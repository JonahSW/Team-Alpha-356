function [t_wait] = Wait(a_1,a_2)
%% Note
% a_1 is the distance from departure planet to sun (AU);
% a_2 is the distance from target planet to sun (AU);

T_1 = 1; % Orbit Period of Departure Planet (years)
T_2 = sqrt((a_2/a_1)^3); % Orbit Period of Target Planet
S = 1/((1/T_1)-(1/T_2)); % Synodic Period
alpha_pi = (((a_1+a_2)/(2*a_1))^(3/2))-1; % time used for wait time in years
y_w = 1:1:5; %iterating years to wait
t_wait = 0;
i = 1;
while t_wait <= 0
    t_wait = (y_w(i) - alpha_pi)*S; % wait time on ceres (years)
    i = i+1;
end

end

