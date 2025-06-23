function V=trid(A,B,C,F,nx)
persistent alpha1 beta1
if isempty(alpha1)
    alpha1=zeros(1,nx);
    beta1=zeros(1,nx);
end
%Forward Thomas path
alpha1(1)=B(1)/C(1);
beta1(1)= F(1)/C(1);
for i=1:nx-1
    zn=1./(C(i)-A(i)*alpha1(i));
    alpha1(i+1)=B(i)*zn;
    beta1(i+1)=(F(i)+A(i)*beta1(i))*zn;
end
%Backward Thomas path
V(nx)=(F(nx)+A(nx)*beta1(nx))/(C(nx)-A(nx)*alpha1(nx));

for i=nx-1:-1:1
    V(i)=V(i+1)*alpha1(i+1)+beta1(i+1);
end
end

