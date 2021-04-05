function [SurTemp,Epsilon,K_MLI] = SurfaceTemp(T,R,D)
%% Note
% T is the surface temperature of the star or the effective tempearture of
% the planet (K);
% R is the radius of the star of the planet (km);
% D is the distance between the Spacecrft and the star/planet (km)
%% Unit Conversion
R = R.*1000; % km -> m
D = D.*1000; % km -> m

%% Constant
MLI_cond = 1E-5; % MLI conductivity (W/m/K)

%% Calculation
SurTemp = T.*sqrt(R./(2.*D));
SurTemp = sum(SurTemp,1);
% toDelete = SurTemp > 1E3;
% SurTemp(:,toDelete) = [];

Epsilon = (6.8E-4).*(SurTemp.^(0.67)); 
Epsilon = sum(Epsilon,1);

% assumed to be the same for all number of layers
K_MLI = MLI_cond.*SurTemp;
K_MLI = sum(K_MLI,1);

end

