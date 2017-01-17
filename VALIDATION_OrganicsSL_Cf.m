% In this code I am assuming that there is no accummulation of organic
% matter beyond the pivot. 
clc;clear all;
Tmax=10;dt=0.001;t=0:dt:Tmax;n=length(t);
%Parameter values
R=0.5;AF=0.5;AS=0.0;
%Positive sea-level fall
ZAratio=0:0.01:3;
% ZAratio=0.1;
zz=length(ZAratio);Cf=zeros(1,zz);
for j=1:zz 
V=0;qorg=0;Lold=0;m=n;Foreset=0;Org=0;Zdot=-ZAratio(j)*AF;
Salt=0;Fresh=0;L=zeros(1,n);U=zeros(1,n);organics=zeros(1,n);
check=zeros(1,n);check1=zeros(1,n);
for i=1:n
Vdot=(1+qorg);V=V+Vdot*dt;
check1(i)=V;
Z=1+Zdot*t(i);
s=(2*V*(1-R))^0.5;r=R*(2*V/(1-R))^0.5;
L(i)=s+Z;
dL=L(i)-Lold;
U(i)=(Z-R*L(i))/(1-R);
if i>2&&dL<0; m=i;LL=L(i);UU=U(i); break;end
Foreset=Foreset+dL*(L(i)-Z);
Salt=Salt+(dL*R-dt*Zdot)*s;
Fresh=Fresh+(dL*R-dt*Zdot)*r;
% check(i)=Foreset+Salt+Fresh;
AccF=-Zdot*r;
AccS=-Zdot*s;
qorg=min(AF*r,AccF)+min(AS*s,AccS);if qorg<0;qorg=0;end
Org=Org+qorg*dt;
organics(i)=Org;
Lold=L(i);
end
for i=m:n
if m==n; break; end
Vdot=(1+qorg);V=V+Vdot*dt;
check1(i)=V;
Z=1+Zdot*t(i);
% Ldot=(Vdot*(1-R)/(LL-Z)+Zdot)/R;
Ldot=(Vdot/(LL-UU)+Zdot)/R;
L(i)=LL+dt*Ldot;
LL=L(i);
U(i)=(Z-R*L(i))/(1-R);
UU=U(i);
if Ldot>0; h=i; break;end
s=L(i)-Z;r=Z-(Z-R*L(i))/(1-R);
Salt=Salt+dt*(-Zdot+R*Ldot)*s;
Fresh=Fresh+dt*(-Zdot+R*Ldot)*r;
% check(i)=Salt+Fresh;
AccF=(-Zdot+R*Ldot)*r;
AccS=(-Zdot+R*Ldot)*s;
qorg=min(AF*r,AccF)+min(AS*s,AccS);if qorg<0;qorg=0;end
Org=Org+qorg*dt;
organics(i)=Org;
end
% Cf(j)=Org/(V-Foreset);
Cf(j)=Org/Fresh;
% Cf(j)=Org/Salt;
end
hold on
plot(ZAratio,Cf,'r','linewidth',2)
% plot(t,L,'r','linewidth',2)
% plot(t,U,'g','linewidth',2)
% plot(t,check1,'g','linewidth',2)
% plot(t,check,'b','linewidth',2)
% xlabel('time');
% ylabel('SH trajectory');
% AXIS([0 Zmax 0 1]);
% legend('R_{ab}=0.2','R_{ab}=0.5','R_{ab}=0.8')
