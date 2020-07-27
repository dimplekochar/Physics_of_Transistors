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
phif=(1/k2)*log(NA/ni);
Eg=1.1;
Vfb=-0.2;
syms psis;
Vdsafine=100;
Vdsa=linspace(0, 2.5, Vdsafine);
Vt=Vfb+(2*phif)+((sqrt(2*epsi*q*NA*2*phif))/Cox);
Vg=Vt+(12*k*T/q); %V 
for m=1:Vdsafine
    Vds=Vdsa(m); %V
    psiSS=vpasolve((Vfb-Vg+psis+(epsi*k1/Cox)*(((k2*psis)+((ni^2/NA^2)*exp(q*(psis)/(k*T))))^0.5))==0, psis);
    psiss=double(psiSS);
    psiSD=vpasolve((Vfb-Vg+psis+(epsi*k1/Cox)*(((k2*psis)+((ni^2/NA^2)*exp(q*(psis-Vds)/(k*T))))^0.5))==0, psis);
    psisd=double(psiSD);
    syms psis;
    f= (un*W/L)*((Cox*(Vg-Vfb-psis))-sqrt(2*epsi*q*NA*psis))*(1+((2*k*T*((Cox*Cox*(Vg-Vfb-psis))+(epsi*q*NA)))/(q*(((Cox*(Vg-Vfb-psis))^2)-(2*epsi*q*NA*psis)))));
    I(m)=vpaintegral(f, psis, [psiss psisd]);
end
plot(Vdsa, I);