%parameters
epsilon=8.854*10^-14; %F/cm^2
epsilonsi=11.8;
epsilonox=3.9;
epsi=epsilon*epsilonsi;
epox=epsilon*epsilonox;
k=1.38*10^-23; %J/K Boltzmann Constant
T=300; %Kelvin
NA=10^10; %cm^-3
ni=10^10;
tox=2*10^-7; %cm
Cox=epox/tox;
q=1.6*10^-19;
k1=q/(k*T); 
k2=ni^2/NA^2;
phif=(1/k1)*log(NA/ni);
Eg=1.1;
Vfb=0;
coeff=sqrt(2*epsi*k*T*NA);

psi=linspace(-0.5,5,200);
Q=coeff.*((((exp(-k1.*psi))+(psi.*k1)-1)+k2.*((exp(k1.*psi))-(psi.*k1)-1)).^0.5);

figure(11); plot(psi,log10(Q));
xlabel("psi");
ylabel("log10(Q)");
grid on;



