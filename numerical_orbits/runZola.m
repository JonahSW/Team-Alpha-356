%% runZola
% sdf33
% run ZolaCharLengths.m
% to iteratively generate a desired mission delta V

%%...initial duration guesses
delVfinal = 20; % km/s, desired mission delta V
delVfinal = delVfinal*1000; % m/s, desired mission delta V
t1 = 500; % days, phase 1 trip time, ACCELERATION

t2_add = 100; % days, additional time on phase 2, COASTING
t2 = t1 + t2_add; % days, phase 2 time step

t3_add = 100; % days, additional time on phase 3, DEACCELERATION
t3 = t2 + t3_add;

Vtol = 0.001; % m/s, velocity tolerance
[delVtot,L,Ttot] = ZolaCharLengths(t1,t2,t3);

if delVtot > delVfinal
    while delVtot > delVfinal
        t1 = t1-1;
        [delVtot,L,Ttot] = ZolaCharLengths(t1,t2,t3);
        fprintf('DeltaV: %.4f km/s\n',delVtot)
        if abs(delVtot-delVfinal)/delVtot < Vtol
            break
        end
    end
else
    while delVtot < delVfinal
        t1 = t1+1;
        [delVtot,L,Ttot] = ZolaCharLengths(t1,t2,t3);
        fprintf('DeltaV: %.4f m/s\n',delVtot)
        if abs(delVtot-delVfinal)*100/delVtot < Vtol
            break
        end
    end
end

fprintf('Mission Totals:\n')
fprintf('   total distance of trip (AU): %.5f\n',L/(1.496e+11))
fprintf('   total duration of trip (days): %.3f\n',Ttot/(3600*24))
fprintf('   total duration of trip (years): %.3f\n',Ttot/(3600*24*365))
fprintf('   total delta V (m/s): %.4f\n',delVtot/1000)
