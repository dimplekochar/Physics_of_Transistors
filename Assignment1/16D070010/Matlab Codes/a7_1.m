%parameters
epsilon=8.854*10^-14; %F/cm^2
epsilonsi=11.8;
epsilonox=3.9;
epsi=epsilon*epsilonsi;
epox=epsilon*epsilonox;
k=1.38*10^-23; %J/K Boltzmann Constant
T=300; %Kelvin
NA=10^17; %cm^-3
ni=1.5*10^10;
tox=5*10^-7; %cm
Cox=epox/tox;
q=1.6*10^-19;
k1=q/(k*T); 
k2=ni^2/NA^2;
phif=(1/k1)*log(NA/ni);
Eg=1.1;
Ei=Eg/2;
Vfb=-(Eg/2+phif);
coeff=sqrt(2*epsi*k*T*NA);
vth=2.6e7; %thermal velocity
sigma=1e-15; %capture cross-section

A1=6e12;
A2=6e12;
A3=2e12;
B1=0.1;
B2=0.97; 
B3=0.5; 
C1=0.2;
C2=0.2;
C3=0.5;

syms eit; %eit=linspace(0,1.1,100);
dit=(A1.*exp(-((eit-B1)./C1).^2))+(A2.*exp(-((eit-B2)./C2).^2))+(A3.*exp(-((eit-B3)./C3).^2));
ditavg=(int(dit, eit, 0, 1.1))/1.1;
dit=ditavg;
tau=(1./(sigma.*ni.*vth)).*exp(-q*(abs(eit-Ei))./(k.*T)); %trap time constant
tauavg=(int(tau, eit, 0, 1.1))/1.1;
tau=tauavg;
w= 10^2 ; %omega

for i = 1:4
    w=w*10;
    Gp=w.*q.*dit.*w.*tau./(1+((w.*tau).^2));
    Cit=q.*dit;
    Rit=tau./Cit;
    psi=linspace(-0.25,1.1,270);

%taking derivative manually and putting the expression
    Cs=-(coeff.*(((exp(-k1.*psi)+psi.*k1-1)+k2.*(exp(k1.*psi)-psi.*k1-1)).^(-0.5)).*((-k1.*exp(-k1.*psi)+k1)+k2.*(k1.*exp(k1.*psi)-k1)))./2;
    Cs(51:270)=-Cs(51:270);
    C1=((1./Cox)+(1./Cs));
    C=1./C1;

%calculating charge
    Q=coeff.*(((exp(-k1.*psi)+psi.*k1-1)+k2.*(exp(k1.*psi)-psi.*k1-1)).^0.5);
    Q(51:270)=-Q(51:270);
    Vg=Vfb+psi-(Q./Cox);

    Cp=Cs+(Cit./(1+((w.*tau).^2)));
    num=(Cox.*((Gp).^2))+(Cp.*((w.*Cox).^2))+(Cox.*((Cp.*w).^2));
    den=((Gp).^2)+((w.^2).*((Cp+Cox).^2));
    Cm=num./den;
%Cmavg=(int(Cm, eit, 0, 1.1)); %/1.1 ?
    Gm=((Cox-Cm).*Gp)./(Cp+Cox);
%Gmavg=(int(Gm, eit, 0, 1.1)); %/1.1 ?

    plot(Vg, Cm);
    hold on;

end;

plot(Vg, C);
legend('w=10^3', 'w=10^4', 'w=10^5', 'w=10^6', 'ideal (no CIT)');
hold off;
xlabel("Vg");
%ylabel("C/Cox"); 
ylabel("C in F");
grid on;