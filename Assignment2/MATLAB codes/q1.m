psisub=V0;
%Evacuum=0
epsilon=8.854*10^-14; %F/cm^2
epsilonsi=11.8;
epsilonox=3.9;
epsi=epsilon*epsilonsi;
epox=epsilon*epsilonox;
tox=5*10^-7; %cm
Cox=epox/tox;
Eg=1.12; %eV
Egox=8;
k=1.38*10^-23; 
T=300;
q=1.6*10^-19;
k11=k*T/q;
NA=10^16;
ni=10^10;
phif=k11*log(NA/ni);
Vfb=-(Eg/2)-phif;
Vt=Vfb+(2*phif)+(sqrt(4*q*NA*epsi*phif)/Cox);
Vg=0.5;
psis=(q*NA*epsi*0.5/Cox^2)*((sqrt(1+(2*Cox^2*(Vg-Vfb)/(q*NA*epsi)))-1)^2);
%gate is 20nm thick, oxide 5nm, substrate 400nm
xgfine=10;
xg=linspace(0, 20*10^-7, xgfine);
xoxfine=20;
xox=linspace(20*10^-7, 25*10^-7, xoxfine);
xsfine=100;
xs(1)=25*10^-7;

if Vg<Vfb %Accumulation

    %poisson for oxide: oxide charge is zero. So, d2psi/dx2=0. So, of the form mx+c.
    Qacc=-Cox*(Vg-Vfb);
    Vox=-Qacc/Cox;
    psiox=(Vox/tox).*(xox-20*10^-7);
    plot(xox, psiox-0.9); %evacuum-eox=0.9
    hold on;
    plot(xox, psiox-0.9-Egox); %oxide bg=8
    borderx1(1)=xox(1);
    borderx1(2)=xox(1);
    borderx2(1)=xox(xoxfine);
    borderx2(2)=xox(xoxfine);
    border1(1)=psiox(1)-0.9;
    border1(2)=psiox(1)-0.9-Egox;
    border2(1)=psiox(xoxfine)-0.9;
    border2(2)=psiox(xoxfine)-0.9-Egox;
    plot(borderx1, border1);
    plot(borderx2, border2);
    %poisson for gate -- we plot fermi level of n+ poly gate
    for i=1:xgfine
        psig(i)=-4.1;
    end
    plot(xg, psig);
    
%poisson for substrate
    %we plot fermi level of substrate
 for i=1:xsfine
        psisubf(i)=-4.15-(Eg/2)-phif;
 end
plot(xscale+xs(1), psisubf-Vg);
plot(xscale+xs(1), psisubf-((Eg/2)-phif)-psisub-Vg); %valence band
plot(xscale+xs(1), psisubf+((Eg/2)+phif)-psisub-Vg); %conduction band

elseif Vg>=Vfb && Vg<Vt %Depletion
 

    %poisson for oxide: oxide charge is zero. So, d2psi/dx2=0. So, of the form mx+c.
    psisdo=(q*NA*epsi/(2*Cox*Cox))*((sqrt(1+(2*Cox*Cox*(Vg-Vfb)/(q*NA*epsi)))-1)^2);
    Qdep=-sqrt(2*q*NA*epsi*psisdo);
    Vox=-Qdep/Cox;
    psiox=(Vox/tox).*(xox-20*10^-7);
    plot(xox, psiox-0.9); %evacuum-eox=0.9
    hold on;
    plot(xox, psiox-0.9-Egox); %oxide bg=8
    borderx1(1)=xox(1);
    borderx1(2)=xox(1);
    borderx2(1)=xox(xoxfine);
    borderx2(2)=xox(xoxfine);
    border1(1)=psiox(1)-0.9;
    border1(2)=psiox(1)-0.9-Egox;
    border2(1)=psiox(xoxfine)-0.9;
    border2(2)=psiox(xoxfine)-0.9-Egox;
    plot(borderx1, border1);
    plot(borderx2, border2);
    
     %poisson for gate -- we plot fermi level of n+ poly gate
   for i=1:xgfine
       psig(i)=-4.1;
    end
    plot(xg, psig);
    %poisson for substrate
    %we plot fermi level of substrate
 for i=1:xsfine
        psisubf(i)=-4.15-(Eg/2)-phif;
 end
plot(xscale+xs(1), psisubf-Vg);
plot(xscale+xs(1), psisubf-((Eg/2)-phif)-psisub-Vg); %valence band
plot(xscale+xs(1), psisubf+((Eg/2)+phif)-psisub-Vg); %conduction band

 
else %Inversion
   

    %poisson for oxide: oxide charge is zero. So, d2psi/dx2=0. So, of the form mx+c.

    Qinv=-Cox*(Vg-Vt);
    Qdep=-sqrt(2*q*NA*epsi*2*phif);
    Vox=-(Qdep+Qinv)/Cox;
    psiox=(Vox/tox).*(xox-20*10^-7);
     plot(xox, psiox-0.9); %evacuum-eox=0.9
    hold on;
    plot(xox, psiox-0.9-Egox); %oxide bg=8
   borderx1(1)=xox(1);
    borderx1(2)=xox(1);
    borderx2(1)=xox(xoxfine);
    borderx2(2)=xox(xoxfine);
    border1(1)=psiox(1)-0.9;
    border1(2)=psiox(1)-0.9-Egox;
    border2(1)=psiox(xoxfine)-0.9;
    border2(2)=psiox(xoxfine)-0.9-Egox;
    plot(borderx1, border1);
    plot(borderx2, border2);
       %poisson for gate -- we plot fermi level of n+ poly gate
   for i=1:xgfine
        psig(i)=-4.1;
    end
    plot(xg, psig);
    
    %poisson for substrate
    %we plot fermi level of substrate
 for i=1:xsfine
        psisubf(i)=-4.15-(Eg/2)-phif;
 end
plot(xscale+xs(1), psisubf-Vg);
plot(xscale+xs(1), psisubf-((Eg/2)-phif)-psisub-Vg); %valence band
plot(xscale+xs(1), psisubf+((Eg/2)+phif)-psisub-Vg); %conduction band


end

        