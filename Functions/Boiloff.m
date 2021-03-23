function [outputArg1,outputArg2] = Boiloff(h_E,h_c)
%% Note
% h_E is the height of LEO
% h_C is the height of LCO

%% Constant
a_Earth = 1; % Distance from Earth to Sun 
a_Ceres = 2.766; % Distance from Ceres to Sun
r_Earth = 6378; % Radius of Earth (km)
r_Ceres = 469.730; % Radius of Ceres (km)
Ba_Earth = 0.306; % Bond albedo of Earth
Ba_Ceres = 0.4; % Bond albedo of Ceres
T_Sun = 5772; % Temperture of the Surface of the Sun (K);
r_Sun = 696340; % Radius of the Sun (km);

%% Unit Conversion
a_LEO = h_E/1000; % m-> km
a_LCO = h_c/1000; % m-> km
a_u = 149597870.691; % AU to km
a_earth = a_Earth*a_u; % AU -> km
a_ceres = a_Ceres*a_u; % AU -> km

%% Variable (Half Way)
[T_Earth] = EffectiveTemp(a_Earth,Ba_Earth); % Effective Temperature of Earth [K]
[T_Ceres] = EffectiveTemp(a_Ceres,Ba_Ceres); % Effective Temperature of Ceres [K]
T_m = [T_Sun; T_Earth; T_Ceres]; % Surface Temperature Matrix [K]

R_m = [r_Sun; r_Earth; r_Ceres]; % Radius Matrix [km]

% Distance from Sun to Spacecraft
D_Sun = [a_earth:1E6:a_ceres]; 
% Distance from Earth to Spacecraft
D_Earth = [a_LEO: ((a_ceres-a_LEO)/(length(D_Sun)-1)) :a_ceres]; 
% Distance from Ceres to Spacecraft
D_Ceres = [a_earth: ((a_LCO-a_earth)/(length(D_Sun)-1)) :a_LCO];
D_m = [D_Sun; D_Earth; D_Ceres]; % Distances Matrix [km]

[SurTemp,Epsilon,K_MLI] = SurfaceTemp(T_m,R_m,D_m);
end

