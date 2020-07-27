%parameters
epsilon=8.854*10^-14; %F/cm^2
epsilonsi=11.8;
epsilonox=3.9;
epsi=epsilon*epsilonsi;
epox=epsilon*epsilonox;
k=1.38*10^-23; %J/K Boltzmann Constant
T=300; %Kelvin
NA=10^17; %cm^-3
ni=10^10;
rox=10^18; %oxide charge density
tox=5*10^-7; %cm
Cox=epox/tox;
Qox=(tox*rox*q)/2;
deltaVg=Qox/Cox;

q=1.6*10^-19;
k1=q/(k*T); 
k2=ni^2/NA^2;
phif=(1/k1)*log(NA/ni);
Eg=1.1;
Vfb=-(Eg/2+phif);
coeff=sqrt(2*epsi*k*T*NA);

psi=linspace(-0.25,1.1,270);

%taking derivative manually and putting the expression
Cs_express=-(coeff.*(((exp(-k1.*psi)+psi.*k1-1)+k2.*(exp(k1.*psi)-psi.*k1-1)).^(-0.5)).*((-k1.*exp(-k1.*psi)+k1)+k2.*(k1.*exp(k1.*psi)-k1)))./2;
Cs_express(51:270)=-Cs_express(51:270);
C1_express=((1./Cox)+(1./Cs_express));
C_express=1./C1_express;


%calculating charge
Q=coeff.*(((exp(-k1.*psi)+psi.*k1-1)+k2.*(exp(k1.*psi)-psi.*k1-1)).^0.5);
Q(51:270)=-Q(51:270);
Vg=Vfb+psi-(Q./Cox);

for i = 2:270
    if (C_express(i) > C_express(i-1))
        Cmin=C_express(i-1);
        C_express(i)=Cmin;
    end
end 
plot(Vg, C_express, Vg-deltaVg, C_express); 
%plot(Vg, C_express./Cox, Vg-deltaVg, C_express./Cox);
xlabel("Vg");
%ylabel("C/Cox"); 
ylabel("C in F");
grid on;
