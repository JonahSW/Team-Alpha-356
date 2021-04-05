%mdp81, jjs280
%03/21/2021
%Code for estimating the tank size and mass for cryogenic propellant storage
%Given a target wet mass and tank diameter, estimates the size and wall thickness of the tanks for a given fluid
%User defines storage pressure
%Tanks modelled as a cylinder with hemispherical endcaps (shape=1) or a sphere (shape=2) 
%Uses refprop to get fluid properties

function [dry_mass, surface_area] = tank_size(wet_mass, diameter, shape, fluid)
    %Hydrogen storage state
    H2_temperature = 
    H2_pressure = 3e5;% [Pa] Guess for hydrogen storage pressure (300kPa)
    H2_density = 
    
    
end


function minimize_surface_area()
    diameter = 1;
    volume = 10;
    height = 0.1:0.1:10;
    SA = (diameter*pi).*height + 4*(diameter/2)^2*pi;
    plot(height, SA)
end