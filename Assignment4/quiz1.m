data;
k1=0.025; %mV kt/q
epsilon=8.89*10^-14; %F/cm
ksi=12;
khfo2=20;
ksio2=4;
ephk=epsilon*khfo2;
epil=epsilon*ksio2;
epsi=epsilon*ksi;
thk=2.5*10^-7; %cm
ni=10^10; %cm^-3
q=1.6*10^-19;
k=1.38*10^-23; %J/K Boltzmann Constant
T=300; %Kelvin

%1. IL thickness
Cox=max(Cpre);
til=epil*((1/Cox)-(thk/ephk));
eot=((til/ksio2)+(thk/khfo2))*ksio2;
%2. Doping
Cmin=min(Cpre);
syms NA;
%phif=-k1*log(NA/ni);
%Wmax=sqrt((2*epsi*2*(phif))/(q*NA));
%Csmin=epsi/Wmax;
fcn = @(NA) (((Cox*Cmin)/(Cox-Cmin))-(epsi/sqrt((2*epsi*2*(k1*log(NA/ni)))/(q*NA))));
NA = fzero(fcn, [10^12, 10^22]);

%3. Cfb, Cmg
LD=sqrt((epsi*k1)/(q*NA));
Csfb=epsi/LD;
Cfb=(Cox*Csfb)/(Cox+Csfb);
phif=k1*log(NA/ni);
Wmg=sqrt((2*epsi*(phif))/(q*NA));
Csmg=epsi/Wmg;
Cmg=(Cox*Csmg)/(Cox+Csmg);

%4. Vfb, Vmg, Vt (ideal)
Vfb_i=-phif; %since midgap workfunction
Vmg_i=Vfb_i+phif+((sqrt(2*q*NA*epsi*phif))/Cox);
Vt_i=Vfb_i+(2*phif)+((sqrt(4*q*NA*epsi*phif))/Cox);

%5. Vfb, Vmg (as fabricated)
t1=0;
t2=0;
for i=1:55
    if ((Cpre(i)<Cfb) && t1==0)
        t1=i;
    end
end
Vfb_af=(Vg(t1)+Vg(t1-1))/2;
for i=1:55
    if ((Cpre(i)<Cmg) && t2==0)
        t2=i;
    end
end
Vmg_af=(Vg(t2)+Vg(t2-1))/2;

%6. Vfb, Vmg (stressed)
t3=0;
t4=0;
for i=1:55
    if ((Cpost(i)<Cfb) && t3==0)
        t3=i;
    end
end
Vfb_s=(Vg(t3)+Vg(t3-1))/2;
for i=1:55
    if ((Cpost(i)<Cmg) && t4==0)
        t4=i;
    end
end
Vmg_s=(Vg(t4)+Vg(t4-1))/2;

%6 fixed charges for as fabricated and stressed device
Qox_af=Cox*(Vmg_af-Vmg_i);
Qox_s=Cox*(Vmg_s-Vmg_i);

%7 dit later

%5, 6 Vt (as fabricated and stressed)
Vt_af=(Cox*Vfb_af+Qox_af+q*4*10^13*phif)/Cox;
Vt_s=(Cox*Vfb_s+Qox_s+q*4*10^13*phif)/Cox;

%8 lfcv, hfcv
k1=0.025; %mV kt/q

coeff=sqrt(2*epsi*k1*q*NA);
k2=ni^2/NA^2;
k3=(1/k1);

psi=linspace(-0.35,1.05,280);
Q=coeff.*((((exp(-k3.*psi))+(psi.*k3)-1)+k2.*((exp(k3.*psi))-(psi.*k3)-1)).^0.5);
Q(71:280)=-Q(71:280);
Vgg=Vfb_i+psi-(Q./Cox);
Cs=-(coeff.*(((exp(-k3.*psi)+psi.*k3-1)+k2.*(exp(k3.*psi)-psi.*k3-1)).^(-0.5)).*((-k3.*exp(-k3.*psi)+k3)+k2.*(k3.*exp(k3.*psi)-k3)))./2;
Cs(71:280)=-Cs(71:280);
C1=((1./Cox)+(1./Cs));
C=1./C1;
t5=0;
for i = 2:280
    if ((C(i) > C(i-1)) && t5==0)
        t5=i;
    end
end
Chf(1:t5-1)=C(1:t5-1);
Chf(t5:280)=C(t5-1);
figure(3);
plot(Vgg, C, Vgg, Chf); 
legend('lfcv', 'hfcv')
xlabel("Vg");
%ylabel("C/Cox"); 
ylabel("C in F");
grid on;

%dit
delpost;
ditpo=-(Cox.*delvgpost-Qox_s)./(q.*(psi_post-phif));
figure(1);
plot(psi_post, ditpo);
delpre;
ditpr=-(Cox.*delvgpre-Qox_af)./(q.*(psi_pre-phif));
figure(2);
plot(psi_pre, ditpr);
