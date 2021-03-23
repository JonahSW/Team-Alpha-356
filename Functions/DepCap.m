function [Target_m0_mf,Target_v,Target_h,Table] = DepCap(r,h,mu,delV,ISP)
%% Note
% r is radius of the departure/capture planet (km);
% h is the height of the departure/capture orbit (km);
% delV is the departure/capture impulsive velocity change (km/s);
% ISP is the total specific impulse of the thrust;
%% Constant
g = 9.807; % gravity m/s^2
EAI = 23.450; % Earth Axis of Inclination (Degrees)
COP = 10.110; % Ceres Orbit Plane (Degrees)
EPOI = 28.573; % Earth Parking Orbit Inclination (Degrees)
DPC = (EAI + COP) - EPOI;
Title = ["Orbital Height(km)"; "Radius + Orbital Height(km)";...
    "Circular Orbit Velocity (km/s)";"Ecentricity of Ellipse (km)";...
    "Semi Minor Axis (km)";"Flight Path Angle at Second Burn (degrees)";...
    "Dep/Cap Velo (km/s)";"Dep/Cap (km/s)";...
    "Velocity Associated with Plane Change (km/s)";"DelV then Plan Change (km/s)";...
    "DelV and Plane Change at the Same Time (km/s)";"Mass Fraction"];
%% Departure
a_0 = h + r; % Departure Orbit = Radius + Height (km)
V_cir = sqrt(mu./a_0); % Circular Orbit Velo (km/s)
a = mu/(delV^2); % Semi Major Axis of Ellipse
e = a_0./a+1; % Ecentricity of Ellipse
b = a.*sqrt(e.^2-1); % Semi Minor Axis
psi = atand(b./a); % Flight Path Angle at 2nd burn (degrees)
V_hyp = (((e+1)./(e-1)).^0.5)*delV; % Departure Velo (km/s)
delV_D = V_hyp-V_cir; % Departure Burn (km/s)
delV_pc = 2.*V_cir.*sind(DPC/2); % Velocity Associated with Plane Change (km/s)
delV_pD = delV_pc + delV_D; % Del V then plane change (km/s)
delV_ppD = sqrt(V_cir.^2+V_hyp.^2-2.*V_cir.*V_hyp.*cosd(DPC));% del V and plane change at the same time (km/s)
m0_mf = exp(1000*delV_D./(g*ISP));% mass fraction

Table = [Title [h;a_0;V_cir;e;b;psi;V_hyp;delV_D;delV_pc;delV_pD;delV_ppD;m0_mf]];
Target_m0_mf = min(m0_mf);
[i] = find(m0_mf == Target_m0_mf);
Target_v = delV_ppD(i);
Target_h = h(i);
end