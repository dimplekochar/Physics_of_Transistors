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
tox=5*10^-7; %cm
Cox=epox/tox;
q=1.6*10^-19;
k1=q/(k*T); 
k2=ni^2/NA^2;
phif=(1/k1)*log(NA/ni);
Eg=1.1;
Vfb=-(Eg/2+phif);

coeff=sqrt(2*epsi*k*T*NA);
psi=linspace(-0.6,1.4,200);
Q=coeff.*(((exp(-k1.*psi)+psi.*k1-1)+k2.*(exp(k1.*psi)-psi.*k1-1)).^0.5);

Vg(1:60)=Vfb+psi(1:60)-(Q(1:60)./Cox); %till psi-s is negative, accumulation, Q is positive
Vg(61:200)=Vfb+psi(61:200)+(Q(61:200)./Cox); %psi-s is positive, Q is negative (Q is substrate charge)

figure(21); 
%plot(psi(1:35),Vg(1:35)); %for more negative psi-s we see exponential increase in negative Vg
plot(psi(35:170),Vg(35:170)); %normal range of psi-s
%plot(psi(160:200),Vg(160:200));  %for more positive psi-s we see exponential increase in Vg
%plot(psi, Vg); %entire plot
xlabel("psi-s");
ylabel("Vg");
grid on;


