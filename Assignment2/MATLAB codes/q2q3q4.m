% cmos.m version 1.3
%
% 1-dimensions selfconsistent Poisson equations, classical solution
% for (100) NMOS or PMOS inversion layer
%
% Originated by Hiroshi Fujioka at U.C. Berkeley, Jun 29 1996
% Revised by S. Kamohara        at U.C. Berkeley, Oct 05 1996
% Revised by Ya-Chin King       at U.C. Berkeley, Feb 11 1997
% Last Revised by Kevin Yang    at U.C. Berkeley, Jun 01 1999
% Posted to http://www-device.eecs.berkeley.edu/~kjyang/qmcv/
%
% Units for Poisson are cgs
% Units for Shrodinger are Hartree (e=m0=h-=1, 1Hr=27.2118eV 1a.u.=0.526177A)
%
% TO RUN THIS PROGRAM PLEASE INCLUDE TWO FUNCTIONS IN THE SAME
% DIRECTORY: shelec.m ,shhole.m
%
%****************initial definition!*********************
%
clear;             		% Clear memory and print header
cflag=0;           		% 0/1 quantum/classical
imax=5000;			% maximum number of newton iterations
%
%****************INPUT PARAMETERS!!!*********************
%
fid1=fopen('cvdata','w');	% Initialize output file; first arg is filename
T = 300;            		% temperature (K)
Ns =1E16;    			% Substrate doping concentration (cm^-3)
				% (-/+ for n/p-type);
Tox=5E-7;         		% Oxide thickness (cm)
pd = 1;            		% consider poly depletion? 1: yes 0: no
Npoly=-5e19;        		% Poly doping (- means n-type; + means p-type)
Vfb=99;				% 99 means autocalculate, else uses value
				% autocalculate uses Npoly and Ns (see below)
				% For metal gate, please input Vfb here
Dvc=  1e-4;       		% small signal voltage for capacitance
np=sign(Ns);			% (do not edit this line) n-sub or p-sub
Vstart=(-1)*np*0.3;       	% Silicon Voltage: start (accumulation)
Vend=np*1.4;         		% Silicon Voltage: end   (inversion)
Nvstep=20;         		% number of voltage steps 
%
cimp0=0;			% Following parameters for implant
%cimp0=3.0e17;			% Gaussian implant: peak conc. (+: p-type)
%ximp0=0.02e-4;			% Gaussian implant: Rp
%simp0=0.05e-4;			% Gaussian implant: dRp (straggle)
%
iii=1;
jjj=1;
%*********DEFINE THE RANGE FOR CALCULATION******************
%
if Nvstep>1;
  Vs0=linspace(Vstart,Vend,Nvstep);
else
  Vs0(1)=Vstart;
end
%
%*********calculate Vfb*****************
%
Ns=abs(Ns);
if Vfb==99
  if np*Npoly<0;
    Vfb=-np*0.026*log(abs(Npoly*Ns)/2.1045e20); 	% n+/nmos or p+/pmos
  else
    Vfb=np*0.026*log(abs(Npoly/Ns)); 	%  n+/nsub or p+/psub
  end
end
%
%****************CONSTANTS**********************
%
k = 8.61735E-5;    % Boltzmann constant in eV/K
beta=1/k/T;        % kT/q (1/eV) 
hb= 6.58215e-16;   % Plank's constant in eV-s
eps0 = 8.86E-14;   % Permittivity of free space (F/cm)
eps1 = 11.7;       % Relative permittivity of Si
eps1ox = 3.9;      % Relative permittivity of Si
q0 = 1.602E-19;    % electron charge (C)
au = 0.5262E-8;    % atomic unit in cm
Eg = 1.12;         % Bandgap of Si in eV
Nv = 1.02E19;      % effective density of states of valence band
Nc = 2.8E19;       % effective density of states of conduction band 
ni = 1.45E10;      % density of intrinsic carrier 
pai= 3.14159;      
m0 = 9.109534E-31; % electron mass in Kg
m1 = 0.916;        % Electron effective mass in z direction for lower(3,6)
m2 = 0.19;         % Electron effective mass in z direction for higher(1,4,2,5)
md1 = 0.19;        % Density of state mass low energy valley (3,6)
md2 = 0.417;       % Density of state mass high energy valley (1,4,2,5)
mc1 = 0.19;        % Conductivity     mass low energy valley (3,6)
mc2 = 0.315;       % Conductivity     mass high energy valley (1,4,2,5)
%
%****************NUMERICAL CONSTANTS******************
%
aregion=2000e-8; % size of classically allowed region (cm)
fregion=100e-8;	% size of classically forbidden region (cm)
N1 =50;		%number of mesh points in the fregion for poisson solution
N2 = 50;	%number of mesh points in the aregion for poisson solution
NS =50;		%number of mesh points in the aregion for Shrodinger eqn.
N=N1+N2;	% total number of mesh points for the poisson solution
dx0= fregion/real(N1);   		% mesh size in fregion
dx1= (aregion-fregion)/real(N2-1);	% mesh size in aregion
%
%**************** MESH GENERATION    **************** 
%
dx=zeros(N,1);			% initialize spacing array
xscale=zeros(N,1);		% initialize space array
dx(1:N1)=ones(N1,1).*dx0;    	% fregion spacing
xscale(2:N1+1)=dx0:dx0:fregion;	% fregion space array
dx(N1+1:N)=ones(N2,1).*dx1;  	% aregion spacing
xscale(N1+1:N)=fregion:dx1:aregion;		% aregion space array
%
%**************** IMPURITY PROFILE SETUP ************************        
%
Nad= Ns*ones(N,1);  		% constant substrate doping
%
if cimp0~=0;			% Gaussian implant
  Nad=Nad+np*cimp0.*exp(-(xscale-ximp0).^2./(simp0^2));
end
%
%****************CALCULATED PARAMETERS**************** 
%
if sign(np)==1
  Nep0=0.5*(Ns+sqrt(Ns*Ns+4*ni*ni));    % hole density at equilibrium
  Nen0=ni*ni/Nep0;                      % electron density at equilibrium
else
  Nen0=0.5*(Ns+sqrt(Ns*Ns+4*ni*ni)); % electron density at equilibrium
  Nep0=ni*ni/Nen0;                      % hole density at equilibrium
end
Ef=+k*T*log(Nep0/ni);                 % definition of Fermi
%
fprintf('Matrix for Poisson is %g by %g square \n',N,N);
Nen= zeros(N,1);  
Nep= zeros(N,1);  
Ne = zeros(N,1);  
Rho = zeros(N,1); 
V0 = zeros(N,1); 
VS = zeros(N,1); 
A = zeros(N,N);  % Matrix for 2nd differential operator 
A(1,1)=1/dx0^2;	 %bondary condition Vsurface=Vs
for j=2:N-1
    avgdx=(dx(j-1)+dx(j))/2;
    A(j,j-1) = 1/dx(j-1)/avgdx;
    A(j,j) = -(1/avgdx)*(1/dx(j-1)+1/dx(j));
    A(j,j+1) = 1/dx(j)/avgdx;
end;
A(N,N)=1/dx(N-1)^2;     
%
%************** BEGIN CALCULATIONS *************************
%
fprintf(fid1,'Vs          Vg          Cac         Lndc        Lnac        Lpdc        Lpac        Es          Qtot        Qinv        \n');
%
%**************LOOP FOR EACH VOLTAGE STEP*************************
%
for ivl=1:Nvstep
  Vss(1)=Vs0(ivl)-Dvc; % for capacitance calculation purpose (dV)
  Vss(2)=Vs0(ivl);
%
%**************LOOP FOR CAPACITANCE CALCULATION*******************
%
  for icc=1:2
    Vs=Vss(icc);
%
%**************POISSON EQUATION SETUP*****************************
%
    Rho(1)=Vs/dx0^2;  %bondary condition on the surface	 
    Rho(N)=0;  %bondary condition V(N)=0
    V0=A\Rho; %inital guess
    deltaNe = zeros(N,1);  
    deltaRho = zeros(N,1); 
    R = zeros(N,1); 
%
%***************** NEWTON LOOP ********************************** 
%
    for i=1:imax;  
%
      % Hole
      if cflag>0.5     		% cflag=1, Classical case
        Nep=+Nep0*exp(-beta*V0); 
      else 			% Quantum   
        xstart=0.0;  		% define quantization region 
        xend=fregion;
        VSH=+V0+Eg/2.0;       	% shift potential
%
	% SOLVING SHRODINGER EQUATIONS 
        [xscaleO,E1h,E2h,E3h,Y1h,Y2h,Y3h,YY1h,YY2h,YY3h,R0]=shhole(xscale,VSH,xstart,xend,NS,Ef,T);
        for j=1:N
          if xscale(j)>=xstart & xscale(j)<=xend;
            Nep(j) = interp1(xscaleO,R0,xscale(j));
          else
            Nep(j)=+Nep0*exp(-beta*V0(j));
          end    			% scale if
        end     		 	% scale for 
      end				% quantum if
%
      % Electron
      if cflag>0.5		% cflag=1, Classical case
        Nen=+ni^2./Nep;
      else			% Quantum
        xstart=0.0;
        xend=fregion;
        VSH=-V0+Eg/2.0;
        [xscaleO,E1e,E2e,Y1e,Y2e,YY1e,YY2e,R0]=shelec(xscale,VSH,xstart,xend,NS,Ef,T);
        for j=1:N
          if xscale(j)>=xstart & xscale(j)<=xend;
            Nen(j) = interp1(xscaleO,R0,xscale(j));
          else
            Nen(j)=+ni^2/Nep0*exp(beta*V0(j));
          end			% scale if
        end			% scale for
      end				% quantum if
%
      Ne=-sign(np)*Nad-Nen+Nep; 	% Net charge density
      deltaNe=-beta*Nep-beta*Nen; 	% gradient for NR method

      Rho=-q0*Ne/eps0/eps1;
      deltaRho=-q0*deltaNe/eps0/eps1;
%
      %Set up the Newton Raphson matrix
      NR=A;
      for j=2:N-1
        NR(j,j)=NR(j,j)-deltaRho(j);
      end 
      R=-A*V0+Rho;	
%
      R(1)=0;
      R(N)=0;
      deltaV0=NR\R;  % Potential in eV
      V0=V0+deltaV0; % Update the V0
%
      % convergence test
      dsort=max(abs(deltaV0));
      if(dsort<1.0e-12)
        break
      end
%
      ivl,i,dsort	              	% print status
    end 				% Newton loop
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% store charge density profile
    if icc<1.5;
      Nen1=Nen;
      Nep1=Nep;
      Ne1=Ne;
    else
      Nen2=Nen;
      Nep2=Nep;
      Ne2=Ne;
    end
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
    Qt(icc)=q0*trapz(xscale,Ne);
    if sign(np)==1
      Qinv(icc)=q0*trapz(xscale,Nen);
    else
      Qinv(icc)=q0*trapz(xscale,Nep);
    end
%
%+++++++++++++++Volatge estimation++++++++++++++++++++++++++++++
%
    Fg0(icc)=(V0(1)-V0(2))/dx(1);                    	% surface electric field
    Vox(icc)=Fg0(icc)*Tox*eps1/eps1ox;               	% potential in oxide
    qtot(icc)=eps0*eps1*Fg0(icc);   	% total charge density (#/cm2)(elec+ion)
    if (Vs*sign(Npoly))<=0; 
      Vpoly(icc)=-0.5*Qt(icc)^2/(eps0*eps1*Npoly*q0);	% poly depetion
    else
      Vpoly(icc)=0.0;					% poly acc
    end
    Vg(icc)=V0(1)+Vox(icc)+Vpoly(icc)*pd+Vfb;         	% gate voltage
    
%
%++++++++++++++++DC Centroid++++++++++++++++++++++++++++++++++++
%
    Lndc1=trapz(xscale,xscale.*(Nen+Nad*((np-1)/2)));   % subtract for pmos 
    Lndc2=trapz(xscale,(Nen+Nad*((np-1)/2)));
    Lndc(icc)=Lndc1/Lndc2;
    Lpdc1=trapz(xscale,xscale.*(Nep-Nad*((np+1)/2)));   % subtract for nmos
    Lpdc2=trapz(xscale,(Nep-Nad*((np+1)/2)));
    Lpdc(icc)=Lpdc1/Lpdc2;
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
  end			% capacitance loop
Vgg(jjj)=Vg(1);
    jjj=jjj+1;
%++++++++++++++++Capacitance++++++++++++++++++++++++++++++++++++
%
  Cac=abs(Qt(2)-Qt(1))/(Vg(2)-Vg(1)) ;
  Cpoly(iii)=Cac;
  iii=iii+1;
%
%++++++++++++++++Ac Centroid++++++++++++++++++++++++++++++++++++
%
  Lnac1=trapz(xscale,xscale.*(-Ne1+Ne2));
  Lnac2=trapz(xscale,Nen1-Nen2);
  Lnac=Lnac1/Lnac2;
  Lpac1=trapz(xscale,xscale.*(Ne1-Ne2));
  Lpac2=trapz(xscale,Nep1-Nep2);
  Lpac=Lpac1/Lpac2;
%
%+++++++++++++++Out put+++++++++++++++++++++++++++++++++++++++++
%
fprintf(fid1,'%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n',Vs,Vg(1),Cac,Lndc(1),Lnac,Lpdc(1), Lpac,Fg0(1),qtot(1),Qinv(1));
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
end			% bias loop
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
plot(Vgg, Cpoly)
fclose(fid1);
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
