%% runZola
% sdf33
% run ZolaCharLengths.m
% to iteratively generate a total duration based on desired delta V
clear; clc; close all;

clear; close all; clc;
disp('low_thrust_modified');
inp=input('Inbound (in) or Outbound (out) -> ','s');
if (strcmp(inp,'out')) || (strcmp(inp,'OUT')) || (strcmp(inp,'Out')) || (strcmp(inp,'o')) || (strcmp(inp,'O'))

%% Heliocentric Outbound
%%...SELECT THESE VALUES...%%
md = 428155.6845; % kg, initial wet mass (Mass after kick stage burn)
mp = 147000; % kg, propellant mass
delVfinal = 1.3646e4; % m/s, desired mission outbound delta V
Vtol = 1e-3; % m/s, velocity tolerance
t1 = 500; % days, phase 1 trip time, ACCELERATION
t2_add = 50; % days, additional time on phase 2, COASTING
t3_add = 80; % days, additional time on phase 3, DEACCELERATION
%%...SELECT THESE VALUES...%%

%%...NOTE ensure that total delta V > desired delta V !!

t2 = t1 + t2_add; % days, phase 2 time step
t3 = t2 + t3_add; % days, phase 3 time step

[delVtot, L, Ttot] = ZolaCharLengths(md,mp,t1,t2,t3);

if delVtot > delVfinal
    while delVtot > delVfinal
        t1 = t1-1;
        t2 = t1 + t2_add;
        t3 = t2 + t3_add;
        [delVtot, L, Ttot] = ZolaCharLengths(md,mp,t1,t2,t3);
        fprintf('DeltaV: %.4f km/s\n',delVtot)
        if abs(delVtot-delVfinal)/delVtot < Vtol
            break
        end
    end
end

fprintf('\n')
fprintf('Outbound Heliocentric Mission Totals:\n')
fprintf('   Total distance of trip (AU): %.5f\n',L/(1.496e+11))
fprintf('   Total duration of trip (days): %.3f\n',Ttot/(3600*24))
fprintf('   Accleration duration of trip (days): %.3f\n',t1)
fprintf('   Decceleration duration of trip (days): %.3f\n',(Ttot)/(3600*24)-t2)
fprintf('   Total duration of trip (years): %.3f\n',Ttot/(3600*24*365))
fprintf('   Total delta V (m/s): %.4f\n',delVtot/1000)


elseif (strcmp(inp,'in')) || (strcmp(inp,'IN')) || (strcmp(inp,'In')) || (strcmp(inp,'i')) || (strcmp(inp,'I')) %From input if statement
%% Heliocentric Inbound

%%...SELECT THESE VALUES...%%
md = 319320; % kg, initial wet mass at Ceres departure
mp = 319320-266852; % kg, propellant mass
delVfinal = 1.3646e4; % m/s, desired mission inbound delta V
Vtol = 1e-3; % m/s, velocity tolerance
t1 = 500; % days, phase 1 trip time, ACCELERATION
t2_add = 50; % days, additional time on phase 2, COASTING
t3_add = 80; % days, additional time on phase 3, DEACCELERATION
%%...SELECT THESE VALUES...%%

%%...NOTE ensure that total delta V > desired delta V !!

t2 = t1 + t2_add; % days, phase 2 time step
t3 = t2 + t3_add; % days, phase 3 time step

[delVtot, L, Ttot] = ZolaCharLengths(md,mp,t1,t2,t3);

if delVtot > delVfinal
    while delVtot > delVfinal
        t1 = t1-1;
        t2 = t1 + t2_add;
        t3 = t2 + t3_add;
        [delVtot, L, Ttot] = ZolaCharLengths(md,mp,t1,t2,t3);
        fprintf('DeltaV: %.4f km/s\n',delVtot)
        if abs(delVtot-delVfinal)/delVtot < Vtol
            break
        end
    end
end

fprintf('\n')
fprintf('Inbound Heliocentric Mission Totals:\n')
fprintf('   Total distance of trip (AU): %.5f\n',L/(1.496e+11))
fprintf('   Total duration of trip (days): %.3f\n',Ttot/(3600*24))
fprintf('   Accleration duration of trip (days): %.3f\n',t1)
fprintf('   Decceleration duration of trip (days): %.3f\n',(Ttot)/(3600*24)-t2)
fprintf('   Total duration of trip (years): %.3f\n',Ttot/(3600*24*365))
fprintf('   Total delta V (m/s): %.4f\n',delVtot/1000)


end%End input if statement