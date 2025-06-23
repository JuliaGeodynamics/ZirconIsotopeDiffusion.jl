function Th=TemperatureHistory(Td,T0,n,amp,t)
tfin=t(end);
Th=max(T0,Td-(Td-T0)*t/tfin-amp*sin(n*2*pi*t/tfin));
% Th = awgn(Tha,0.1,0.1);
%XH2O=3.*(1+0.25*sin(30.*2*pi*tt/tfin));

% 
% %Dyke
% h=20; %half thickness of the dike
% x=0; % coordinate for T evolution
% k=0.6e-6; %temeparature conduction coeff in m2/s
% Th=(Td - T0) * (erf((h - x) * (k * tt)^(-0.5)/2)/2 + erf((h + x) * (k * tt)^(-0.5)/2)/2) + T0;
end