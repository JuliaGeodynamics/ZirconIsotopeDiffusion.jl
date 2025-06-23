clear all
nt=1000; nx=1000; %number of points in time and space
figure(1); clf;  figure(2);clf; pause(1e-10);
xi=1-cos(pi*linspace(0,1,nx)/2).^2; A=0*xi; B=A; C=A; F=A;
V=0; W=0; tyear=365*3600*24;
tfin=5000*tyear; dt=tfin/(nt-1); t=0:dt:tfin; %time array
Td=840+273; %T of the dike
T0=720+273; %T of country rocks
c=0*A+ZrSaturation(Td); Czirc=490000; kpl=0.2; XH2O=2;
Th=TemperatureHistory(Td,T0,10,3,t);
figure(4); plot(t/tyear,Th-273); xlabel ('time'), ylabel('Temperature, ^oC'); pause(1e-5);
%Trace element initial concentrations
CH0(1:nx,1)=4.7;%Hf  %FCT from Bachmann etal JPet 2002.
CH0(1:nx,1)=6;%Hf  %FCT from Bachmann etal JPet 2002.
CH0(1:nx,2)=21.9; %Y
CH0(1:nx,3)=3.9;%U
CH0(1:nx,4)=12.4; %Th
CH0(1:nx,5)=4.;%Sm;%P
CH0(1:nx,6)=4.;%Dy
CH0(1:nx,7)=4;%Yb
% CH0(1:nx,8)=1000;%P
ntr=size(CH0,2);
CH=zeros(1,nx);
CHS=zeros(nt,ntr);
CHS(1,:)=CH0(1,:).*Kdd(Th(1));

Rc=0.1; s=1e-4; %cell and Zircon radius in cm;
Mas0=4*pi*s^3*Czirc/3+4*pi*c(1)*(Rc^3-s^3)/3; Mpl=0;
Dscale=DiffusionCoefficient((Th(1)+Th(nt))/2,3.0);
tscale=Rc^2/Dscale; dt=dt/tscale;
t=t/tscale; s=s/Rc; R=1;
for j=2:nt %loop in time
    c0=c; s0=s; R0=R; errC=1e10; itt=0;
    Di=DiffusionCoefficient(Th(j),XH2O)/Dscale; cm=ZrSaturation(Th(j));
    dXpdT=(Xplag(Th(j)+1e-5)-Xplag(Th(j)-1e-5))/2e-5;  dTdt=(Th(j)-Th(j-1))/dt;
    while errC>1e-7 && itt<1000 %Zircon growth
        ci=c;
        V=-Di*(ci(2)-cm)/(xi(2)-xi(1))/(R-s)/(cm-Czirc); %Positive when zircon is growing
        W=-dXpdT*dTdt/3/R^2; %Negative when major phase growth.
        R=R0+W*dt; s=s0+V*dt;  H=R-s; ktr=0;
        [A,B,C,F]=coefs(A,B,C,F,c0,xi,H,V,W,R,s,cm,kpl,ktr,nx,Di,dt,0);
        c=trid(A,B,C,F,nx); errC=max(abs(c-ci)); itt=itt+1;
    end
    kt=Kdd(Th(j)); Dif=TraceDiff(Th(j),XH2O)/Dscale;
    for i=1:ntr
        kplt=0;  ktr=kt(i);  Di=Dif(i);
        [A,B,C,F]=coefs(A,B,C,F,CH0(:,i),xi,H,V,W,R,s,cm,kplt,ktr,nx,Di,dt,i);
        CH=trid(A,B,C,F,nx); CHS(j,i)=ktr*CH(1); CH0(:,i)=CH;
    end
    Mpl=Mpl-4*pi*R^2*W*dt*kpl*c(nx); Xs(j)=s*Rc*1e4; if V<0, disp('dissolution'); end
    if mod(j,nt/100)==0
        figure(1);r=(s+(R-s)*xi)*Rc; V0=Rc*1e4/tscale*tyear;
        Mass=4*pi*(s*Rc)^3*Czirc/3+4*pi*trapz(r,r.^2.*c)+Mpl*Rc^3;
        nexttile(1); plot(r*1e4,c); xlabel('r, \mum'); ylabel('C, ppm'); hold on
        nexttile(2); plot(t(j)*tscale/tyear,(Mass-Mas0)/Mas0,'+r'); xlabel('time, years'); ylabel('\deltaM'); hold on
        nexttile(3); plot(t(j)*tscale/tyear,V*V0,'+r',t(j)*tscale/tyear,-W*V0,'ob'); xlabel('time, years'); ylabel('V,-W, \mum/year');hold on
        nexttile(4); plot(t(j)*tscale/tyear,s*Rc*1e4,'+r'); xlabel('time, years'); ylabel('Zr rad, \mum'); hold on
        pause(1e-11)
    end

end
figure(2)
CHmax=max(CHS(:,1:ntr));
plot(Xs(2:nt),CHS(2:nt,1:ntr)./CHmax,LineWidth=2);
legend(["Hf","Y","U","Th","Sm","Dy","Yb"])
