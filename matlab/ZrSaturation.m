function Csat=ZrSaturation(T)
%Csat=4.414e7/exp(13352/T)/2; %Watson 96, Eq 1, in ppm Zr (divide by 1),or mol Zr (divide by 2) 
Mfactor = 0.0000048*(T)^2-0.0083626*(T)+4.8484463; % empirical relations from magma
%differentiation calc (file  M_factorsforOleg.xlsx
% Mfactor=1.32;
Csat=490000/exp(10108/T-1.16*(Mfactor-1)-1.48); % Boehnkeetal2013ChemGeol351,324 assuming 490000ppm in Zircon
%Csat=490000/(exp(12900/T-0.85*(Mfactor-1)-3.80));
end