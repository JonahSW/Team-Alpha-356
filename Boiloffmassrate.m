function [R_S_P, BOM_Rate] = Boiloffmassrate(BOM_important, time_important, R_s_p)

a_u = 1.5E+8; % AU to Km
R_S_P = (a_u*1000)*(R_s_p((30:length(R_s_p))));
BOM_Rate= real(BOM_important./time_important);

end 