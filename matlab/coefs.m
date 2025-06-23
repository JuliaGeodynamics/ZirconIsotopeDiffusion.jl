function [A,B,C,F]=coefs(A,B,C,F,c0,xi,H,V,W,R,s,cm,kpl,ktr,nx,Di,dt,itrace)
if itrace ==0 %boundary conditions for zirconium
    A(1)=0;B(1)=0;C(1)=1;F(1)=cm;
else %for trace element
    C(1)=Di-V*(xi(2)-xi(1))*(R-s)*(1-ktr); A(1)=0; B(1)=Di; F(1)=0;
end

for i=2:nx-1
    xi12=(xi(i-1)+xi(i))/2; h12=xi(i)-xi(i-1);
    xi23=(xi(i+1)+xi(i))/2; h23=xi(i+1)-xi(i);
    t3 = (H * (V - W));
    t4 = (xi12 ^ 2);
    t10 = (2 * V - W) * s - V * R;
    t12 = V * s;
    t25 = (xi23 ^ 2);
    t40 = H ^ 2;
    A(i) = -(1 / h12 * dt * (h12 * H * (xi12 * t10 + t4 * t3 - Di - t12) + 2 * (xi12 * H + s) * Di)) / 0.2e1;
    B (i)= (1 / h23 * (h23 * H * (xi23 * t10 + t25 * t3 - Di - t12) - 2 * (xi23 * H + s) * Di) * dt) / 0.2e1;
    mt = (((-xi12 + 2 - xi23) * s + R * (xi12 + xi23)) * t40 * (xi12 - xi23)) / 0.2e1;
    C(i)=A(i)+B(i)+mt;
    F(i)=mt*c0(i);
end
C(nx)=Di+W*(xi(nx)-xi(nx-1))*(R-s)*(1-kpl); A(nx)=Di; B(nx)=0;
end