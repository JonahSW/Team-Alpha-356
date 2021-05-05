%% Zola Characteristic Lengths
function [delVtot, L, Ttot] = ZolaCharLengths(md,mp,t1,t2,t3)
% 5.2.21
% sdf33

%...Project Demeter
%...low-thrust trajectory from Earth to Ceres, and return
%...constant thrust
%...note: t_ , time, values start from t = 0

%...general parameters
g = 9.807; % m/s2, Earth gravity

%...Alpha parameters
mo = md + mp; % kg, wet mass
F = 120.1; % N, thrust
ao = F/mo; % m/s2, initial acceleration
Isp = 8180; % sec, specific impulse
vj = Isp*g; % effective velocity

%%...phase 1...accelerating, 0 < t < t1
t1 = t1*24*3600; % secs, time
%a1 = (ao)./(1 - (ao)*t1/vj); % m/s2, acceleration magnitude
V1 = -vj*log(1-(ao)*t1/vj); % m/s, velocity magnitude
L1 = (vj^2/(ao))*((1-(ao)*t1/vj)*log(1-(ao)*t1/vj) + (ao)*t1/vj); % m, characteristic length

% fprintf('Phase 1 Totals:\n')
% fprintf('   total distance of trip (AU): %.5f\n',L1/(1.496e+11))
% fprintf('   total duration of trip (days): %.3f\n',t1/(3600*24))
% fprintf('   total delta V (km/s): %.4f\n',V1/1000)

%%...phase 2...coasting, t1 < t < t2
t2 = t2*24*3600; % secs, time
Vmax = V1; % m/s, maximum velocity during coasting
L2 = Vmax.*(t2-t1); % m, characteristic length

% fprintf('Phase 2 Totals:\n')
% fprintf('   total distance of trip (AU): %.5f\n',L2/(1.496e+11))
% fprintf('   total duration of trip (days): %.3f\n',t2-t1)
% fprintf('   total delta V (km/s): %.4f\n',0)

%%...phase 3...deccelerating, t2 < t < T
T = t3*24*3600; % secs, time
%L3 = vj.^2./ao.*((ao.*t1./vj).^2-(1-ao.*t1./vj).*log(1-ao.*t1./vj) - ao.*t1./vj); % m, characteristic length
L3 = ao*t1^2 - L1;

% fprintf('Phase 3 Totals:\n')
% fprintf('   total distance of trip (AU): %.5f\n',L3/(1.496e+11))
% fprintf('   total duration of trip (days): %.3f\n',t3-t2-t1)
% fprintf('   total delta V (km/s): %.4f\n',V1/1000)

%%...mission totals
L = L1+L2+L3; % m, total distance of trip
Ttot = T; % sec, total duration of trip
delVtot = 2*Vmax; % m/s, total delta V
% fprintf('Mission Totals:\n')
% fprintf('   total distance of trip (AU): %.5f\n',L/(1.496e+11))
% fprintf('   total duration of trip (days): %.3f\n',Ttot/(3600*24))
% fprintf('   total delta V (km/s): %.4f\n',delVtot/1000)

end %...ZolaCharLengths
