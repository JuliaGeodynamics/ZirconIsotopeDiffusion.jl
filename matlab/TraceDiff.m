function Dif=TraceDiff(T,x)
%lnD_Hf=(-8.620340372*T-42705.17449-.318918919*x*T+4049.500765*x)/T;
lnD_Hf=log(DiffusionCoefficient(T,x)/3e4);
lnD_Y=-9.925340370-(35587.24428-2615.214492*x)/T;
P=0.01; %Pressure in GPa
lnD_U=-6.37-2.65*sqrt(x) - (44729-1093*P-8944*sqrt(x))/T;
lnD_Th=-7.02-1.30*sqrt(x) - (44682-1370*P-7281*sqrt(x))/T;
lnD_P=((0.9469e1 * x - 0.1815e3) / (x + 0.2374e1) * ((1000 / T) - 0.54e0) - 0.267e2);
Dif=exp([lnD_Hf,lnD_Y,lnD_U,lnD_Th,lnD_Y,lnD_Y,lnD_Y])*1e4; %,lnD_P
end