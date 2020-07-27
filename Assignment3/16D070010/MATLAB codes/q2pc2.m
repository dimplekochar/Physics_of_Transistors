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
Vdsfine=200;
Vds=linspace(0, 2.5, Vdsfine);
Vg=2*(k*T/q)*log(NA/ni); %V Vg=Vth
Vt=2*(k*T/q)*log(NA/ni);
Vox=((sqrt(2*epsi*q*NA*2*phif))/Cox);
m=1+(sqrt(epsi*q*NA/(4*phif))/Cox);
Isat=(un*Cox*W/L)*((Vg-Vt)^2)/(2*m);
Imax=0;
if Vg>Vt
for i=1:Vdsfine
    I(i)=(un*Cox*W/L)*(((Vg-Vfb-(2*phif)-(Vds(i)/2))*Vds(i))-((2*sqrt(2*epsi*q*NA)/(3*Cox))*(((2*phif+Vds(i))^1.5-(2*phif)^1.5))));
   if I(i)>=Imax
        Imax=I(i);
    else
        I(i)=I(i-1); %I(i)=Isat;
    end
end

else 
    for i=1:Vdsfine
    I(i)=(un*Cox*W*(m-1)/(L*k2*k2))*((exp(k2*(Vg-Vt)/m))*(1-(exp(-k2*Vds(i)))));
    end
    
    
end
plot(Vds, I);
