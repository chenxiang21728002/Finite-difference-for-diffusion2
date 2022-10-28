clear; close all;
Dox0=2.29*10^-5;%Diffusion coefficient in pure water (cm2/s)
Dg0=6.75*10^-6;%Diffusion coefficient in pure water (cm2/s)
Da0=6.75*10^-6;%Diffusion coefficient in pure water (cm2/s)
Din0=6.75*10^-6;%Diffusion coefficient in pure water (cm2/s)
Cenz=500*10^-9;%Concentration of glucose oxidase (mol/cm3)
Vmax0=860*Cenz;%Kinetics constants(mol/(s cm3))
C0ox=0.274*10^-6;%Species concentration in the external medium (mol/cm3)
C0g=5.5*10^-6;
C0a=0;
Kox0=6.992*10^-5;%Kinetics constants(mol/cm3)
Kox=Kox0/C0ox;
Kg0=6.187*10^-7;%Kinetics constants(mol/cm3)
Kg=Kg0/C0ox;
l=0.025;%Thickness of the memberme(cm)
L=l/l;
h=0.01;%spatial step
dt=4*10^-5;%s
N=L/h+1;
x=0:h:L;
Psi=zeros(length(x),1);
Lam=zeros(length(x),1);
PH1=zeros(length(x),1);
C1ox=zeros(length(x),1);
C2ox=zeros(length(x),1);
C1oxt=zeros(length(x),1);
C1g=zeros(length(x),1);
C2g=zeros(length(x),1);
C1gt=zeros(length(x),1);
C1a=zeros(length(x),1);
C2a=zeros(length(x),1);
C1at=zeros(length(x),1);
PH2=zeros(length(x),1);
xi=zeros(length(x),1);
C2a=zeros(length(x),1);
Doxtemp=zeros(length(x),1);
Dgtemp=zeros(length(x),1);
Datemp=zeros(length(x),1);
Dintemp=zeros(length(x),1);
Doxtemp1=zeros(length(x),1);
Dgtemp1=zeros(length(x),1);
Datemp1=zeros(length(x),1);
Dintemp1=zeros(length(x),1);
C1ox(:,1)=C0ox/C0ox;
% C2ox(1)=-0.1;
% C2ox(end)=0;
C1g(:,1)=C0g/C0ox;
% C1g(1,1)=C0g;
% C2g(1)=-0.1;
% C2g(end)=0;
C1a(:,1)=C0a/C0ox;
% C2a(1)=0.1;
% C2a(end)=0;
c10=3.59*10^-9;
c20=1.758*10^-6;
c30=0.01774*10^-3;
c40=0.0972*10^-3;
c50=3.98*10^-11;
c60=0.154*10^-3;
Psi0=0.63;%Initial polymer volume fraction
Sigma=0.93*10^-3;%Molar density of the charged group in the polymer (equiv/cm3)
PKg=7.8;%Acid dissociation constants
PK=6.39;
PH0=7.409;%PH of external solution
A=10^(PH0-PKg);
z=0;

for n=1:N
Doxtemp(n)=Dox0;
Dgtemp(n)=Dg0;
Datemp(n)=Da0;
Dintemp(n)=Din0;
Psi(n)=Psi0;
end

for n=1:N
    
a5=(1-Psi(n))*(c50+c60);
a4=(1-Psi(n))*(c50+c60)*A+Sigma*Psi(n);
a3=-(c10+c40)*(1-Psi(n));
a2=-(2*c20+A*c10+A*c40)*(1-Psi(n));
a1=-(3*c30+2*A*c20)*(1-Psi(n));
a0=-(3*c30*A)*(1-Psi(n));
P=[a5,a4,a3,a2,a1,a0];
Lamt=roots(P);
i=find(Lamt>0 & (imag(Lamt)==0)); 
Lam(n)=Lamt(i);
end

for n=1:N
PH1(n)=PH0-log10(Lam(n));
end

c1=(Lam(1)^-1)*c10;
c2=(Lam(1)^-2)*c20;
c3=(Lam(1)^-3)*c30;
c4=(Lam(1)^-1)*c40;
c5=(Lam(1)^1)*c50;
c6=(Lam(1)^1)*c60;

buffer=c1+c2+c3;
B=10^(PH1(1)-PK);

T=0;

error4=1;
while(error4>10^(-15)  | error5>10^(-15)  | error6>10^(-15)| error7>10^(-15));
    
error1=1;

while(error1>10^(-15)  | error2>10^(-15)  | error3>10^(-15));

    Doxa=mean(Doxtemp);
    Dga=mean(Dgtemp);
    Daa=mean(Datemp);
    T=T+1;
    M=T/dt;
    tox=T*l*l/Doxa;
    tg=T*l*l/Dga;
    ta=T*l*l/Daa;
    Vmaxox=Vmax0/C0ox*tox;
    Vmaxg=Vmax0/C0ox*tg;
    Vmaxa=Vmax0/C0ox*ta;
    j=0;
    while(j<=M);
    
    C1oxt(1)=C0ox/C0ox;
    C1gt(1)=C0g/C0ox;
    C1at(1)=C0a/C0ox;
for i=2:length(x)-1;    
   
    
    C1oxt(i)=C1ox(i)+dt/h/h*(C1ox(i+1)-2*C1ox(i)+C1ox(i-1))-dt/T*Vmaxox*C1ox(i)*C1g(i)/2/(C1ox(i)*(Kg+C1g(i))+C1g(i)*Kox);
    C1gt(i)=C1g(i)+dt/h/h*(C1g(i+1)-2*C1g(i)+C1g(i-1))-dt/T*Vmaxg*C1ox(i)*C1g(i)/(C1ox(i)*(Kg+C1g(i))+C1g(i)*Kox);
    C1at(i)=C1a(i)+dt/h/h*(C1a(i+1)-2*C1a(i)+C1a(i-1))+dt/T*Vmaxa*C1ox(i)*C1g(i)/(C1ox(i)*(Kg+C1g(i))+C1g(i)*Kox);
         
end
    C1oxt(end)=C1oxt(end-1);
    C1gt(end)=C1gt(end-1);
    C1at(end)=C1at(end-1);

C1ox1=C1ox;
C1g1=C1g;
C1a1=C1a;

C1ox=C1oxt;
C1g=C1gt;
C1a=C1at;

j=j+1;
% subplot(3,2,1);plot(x,C1ox);ylabel('C1ox');
%  axis([x(1) x(end) 0 1])
% 
% subplot(3,2,2);plot(x,C1g);ylabel('C1g');
%  axis([x(1) x(end) 0 20])
% 
% subplot(3,2,3);plot(x,C1a);ylabel('C1a');
% axis([x(1) x(end) 0 10])
% getframe;
    end
    
    
error1=max(abs(C1ox-C1ox1));
error2=max(abs(C1g-C1g1));
error3=max(abs(C1a-C1a1));

% subplot(3,2,1);plot(x,C1ox);ylabel('C1ox');axis([x(1) x(end) 0 0.3*10^-6])
% 
% subplot(3,2,2);plot(x,C1g);ylabel('C1g');axis([x(1) x(end) 0 6*10^-6])
% 
% subplot(3,2,3);plot(x,C1a);ylabel('C1a');axis([x(1) x(end) 0 0.3*10^-6])
% getframe;

% subplot(3,2,4);plot(x,PH2);ylabel('PH2');
% subplot(3,2,5);plot(x,Doxtemp/Dox0);ylabel('Doxtemp/Dox0');
% figure(1);plot(x,C1ox);ylabel('C1ox');
% 
%     axis([x(1) x(end) 0 3*10^-7])
%     getframe;

end
% tse=T*L*L/Dse;
% tpe=T*L*L/Dpe;
% x1=0:h*L*10^7:L*10^7;%nm
C1ox=C1ox*C0ox;%mol/cm3
C1g=C1g*C0ox;%mol/cm3
C1a=C1a*C0ox;%mol/cm3



for n=1:N

    t=B-C1a(n)*(1+B)/buffer;
    PH2(n)=PK+log10((B-C1a(n)*(1+B)/buffer)/(1+C1a(n)*(1+B)/buffer));
end

for n=1:N
Lam(n)=10^(PH0-PH2(n));
end

for n=1:N

    D=-c10*Lam(n)^-1-2*c20*Lam(n)^-2-3*c30*Lam(n)^-3-c40*Lam(n)^-1+c50*Lam(n)+c60*Lam(n);
    Psi(n)=D/(D-Sigma/(1+A*Lam(n)^-1));
    
end
CN=14.4;
lc=1.54;%A
rox=1.21;
rg=2.5;
ra=2.5;
rin=16;
Np=100;
for n=1:N

    xi(n)=CN^0.5*lc*Np^0.5*(1/Psi(n))^(1/3);
    Doxtemp1(n)=(1-rox/xi(n))*exp(-1/(1/Psi(n)-1))*Dox0;
    Dgtemp1(n)=(1-rg/xi(n))*exp(-1/(1/Psi(n)-1))*Dg0;
    Datemp1(n)=(1-ra/xi(n))*exp(-1/(1/Psi(n)-1))*Da0;
    Dintemp1(n)=(1-rin/xi(n))*exp(-1/(1/Psi(n)-1))*Din0;
end
 
error4=max(abs(Doxtemp-Doxtemp1));
error5=max(abs(Dgtemp-Dgtemp1));
error6=max(abs(Datemp-Datemp1));
error7=max(abs(Dintemp-Dintemp1));

    Doxtemp=Doxtemp1;
    Dgtemp=Dgtemp1;
    Datemp=Datemp1;
    Dintemp=Dintemp1;
    z=z+1;
%     figure(1);plot(x,C1ox);ylabel('C1ox');

end

subplot(3,2,1);plot(x,C1ox);ylabel('C1ox');
% axis([x(1) x(end) 0 0.3*10^-6])

subplot(3,2,2);plot(x,C1g);ylabel('C1g');
% axis([x(1) x(end) 0 6*10^-6])

subplot(3,2,3);plot(x,C1a);ylabel('C1a');
% axis([x(1) x(end) 0 0.8*10^-6])

subplot(3,2,4);plot(x,PH2);ylabel('PH2');
subplot(3,2,5);plot(x,Doxtemp/Dox0);ylabel('Doxtemp/Dox0');
