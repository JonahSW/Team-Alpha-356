function [Temp_2] = SurfaceTemp(T,R,D)
%% Note
% T is the surface temperature of the star or the effective tempearture of
% the planet (K);
% R is the radius of the star of the planet (km);
% D is the distance between the Spacecrft and the star/planet (km)
%% Unit Conversion
R = R*1000; % km -> m
D = D*1000; % km -> m
%% Calculation
Temp_2 = T*sqrt(R/(2*D));
end

