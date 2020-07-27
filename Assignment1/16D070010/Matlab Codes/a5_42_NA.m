%parameters
epsilon=8.854*10^-14; %F/cm^2
epsilonsi=11.8;
epsilonox=3.9;
epsi=epsilon*epsilonsi;
epox=epsilon*epsilonox;
k=1.38*10^-23; %J/K Boltzmann Constant
T=300; %Kelvin

ni=10^10;
tox=5*10^-7; %cm
Cox=epox/tox;
q=1.6*10^-19;
k1=q/(k*T); 

NA=10^14; %cm^-3
for i=1:3
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

plot(Vg, C_express); 
%plot(Vg, C_numerical./Cox); 
hold on;

%depletion approximation

t=0;
for i=1:270
    if(psi(i)>2*phif) && (t==0)
        t=i-1;
    end
end


depw=sqrt(2.*epsi.*psi(51:t)./(q.*NA));
depww=sqrt(4.*epsi.*phif./(q.*NA));
Cmin=epsi/depww;

%accumulation
Q(1:50)=coeff.*(exp(-k1.*psi(1:50)./2));
C(1:50)=Cox;

%depletion
Q(51:t)=-coeff.*((psi(51:t).*k1).^0.5);
Cd=epsi./depw;
C1=(1./Cox)+(1./Cd);
C(51:t)=(1./C1);

%inversion
Q(t+1:270)=coeff.*(((psi(t+1:270).*k1)+k2.*(exp(k1.*psi(t+1:270)))).^0.5);
C(t+1:270)=Cmin;

Vg(1:50)=Vfb+psi(1:50)-(Q(1:50)./Cox); %till psi-s is negative, accumulation, Q is positive
Vg(51:270)=Vfb+psi(51:270)+(Q(51:270)./Cox); %psi-s is positive, Q is negative (Q is substrate charge)
NA=NA*100;
plot (Vg, C);
end

legend('ideal=1e14', 'dep=1e14', 'ideal=1e16', 'depl=1e16', 'ideal=1e18', 'dep=1e18');
hold off;

xlabel("Vg");
ylabel("C (in F)"); 
grid on;


