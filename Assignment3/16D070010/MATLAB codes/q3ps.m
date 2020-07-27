clear
k=1.38*10^-23; %
T=300; %temperature in kelvin
NA=5*10^15;
epsilon=8.854*10^-14;
ksi=11.8;
kox=3.9;
epsi=ksi*epsilon;
epox=kox*epsilon;
k1=sqrt((2*k*T*NA)/(epsi));
ni=10^10;
q=1.6*10^-19;
tox=5*10^-7;
L= 1;
W=40;
un=800;
Cox=epox/tox;
k2=q/(k*T); 
k3=(k*T)/q;
phif=(1/k2)*log(NA/ni);
Eg=1.1;
Vfb=-0.2;
Vt=Vfb+(2*phif)+((sqrt(2*epsi*q*NA*2*phif))/Cox);
Vg=Vt; %V 
syms psis;
Vdsafine=100;
Vdsa=linspace(0, 3, Vdsafine);
for mm=1:4
for m=1:Vdsafine
Vds=Vdsa(m); %V %fix a vds

vdsfine=200;
V=linspace(0, Vds, vdsfine); %divide vds in vdsfine portions
I(m)=0;

for i=1:vdsfine
    I1(i)=0;
    psiS(i)=vpasolve((Vfb-Vg+psis+((epsi*k1/Cox)*(((k2*psis)+((ni^2/NA^2)*exp(q*(psis-V(i))/(k*T))))^0.5)))==0, psis);
    psi_s(i)=double(psiS(i));
    
    psifine=200;
    psi=linspace(0.01, psi_s(i), psifine);
    for j=1:psifine
        I1(i)=I1(i)+(((q*un*W/L)*((ni*ni/NA)*(exp(q*(psi(j)-V(i))/(k*T))))/((sqrt((2*k*T*NA)/epsi))*(((q*psi(j)/(k*T))+((ni^2/NA^2)*exp(q*(psi(j)-V(i))/(k*T))))^0.5)))*(psi_s(i)/psifine));
    end
    I1(i)=I1(i)*(Vds/vdsfine);
    I(m)=I(m)+I1(i);
end
end
plot(Vdsa, I);
hold on;
Vg=Vg+(5*k3);
end
