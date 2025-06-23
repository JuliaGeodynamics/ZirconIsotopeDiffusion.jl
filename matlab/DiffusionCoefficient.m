function Dif=DiffusionCoefficient(T,x) %X-wt%water
theta=1000/T;
lnD=-(11.4*x+3.13)/(0.84*x+1)-(21.4*x+47)/(1.06*x+1)*theta;%best fit of Zr Diff coefficients (several workers) and WH dependence on XH2O
Dif=exp(lnD)*1e4;% in cm2/s
end
