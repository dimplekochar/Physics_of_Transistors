function [xscale,E1,E2,E3,Y1,Y2,Y3,YY1,YY2,YY3,R]=shhole(xscaleI,VI,xstart,xend,N,Ef,T);
%***********************************************************
% 1-dimensions Shrodinger equations
% for (100) hole ver 1.0.0
% "input"
% xscaleI:cordinate of volatge (cm)
% VI: voltage (V)
% (xstrat,xend): analysis region (cm)
% N: analysis mesh number
% Ef: Fermi level
% T : temperature
% "output"
% xscale : new scale (cm)
% E1-3 (i):energy level of light, heavy, split on i-band
% YY1-3   (j,i):density of light, heavy, split on i-band at j-point
% R : carrirer density (cm-3)
% Written by Hiroshi Fujioka at U.C.Berkeley, Jun 29 1996
% Revised by Shiro Kamohara  at U.C.Berkeley, Dec 06 1996
% Unit for Shrod is Hartree (e=m0=h-=1, 1Hr=27.2118eV 1a.u.=0.526177A)
%****************CONSTANTS*********************************************
k = 8.61735E-5;    % Boltzmann constant in eV/K
hb= 1.05459E-27;   % Plank's constant in ergs
eps0 = 8.86E-14;   % Permittivity of free space (F/cm)
eps1 = 11.7;       % Relative permittivity of Si
eps1ox = 3.9;      % Relative permittivity of Si
q0 = 1.602E-19;    % electron charge (C)
au = 0.5262E-8;    % atomic unit in cm
Hr = 27.212;       % 1 Hartree in eV
Eg = 1.12;         % Bandgap of Si in eV
Nc = 2.8E19;       % effective density of states for conduction band
Nv = 1.02E19;      % effective density of states for valence band
m0 = 9.109534E-28; % electron mass in g
m1 = 0.291;        % Electron mass in z direction for heavy hole
m2 = 0.2;          % Electron mass in z direction for light hole
m3 = 0.29;         % Electron mass in z direction for split off
md1 = 0.645;       % Density of state mass for haevy hole
md2 = 0.251;       % Density of state mass for light hole
md3 = 0.29;        % Density of state mass for sprit-off band
deltaE=0.044;      % spin-orbit sprit energy
g1 = 1;            % degeneracy for heavy hole 
g2 = 1;            % degeneracy for light hole
g3 = 1;            % degeneracy for split off band 
%****************scale set up         ****************
xscale = linspace(xstart,xend,N).';  % New scale cm
dx0=(xend-xstart)/real(N);
dx= dx0/au;             % Mesh separation in a.u.
dd=1/2/(dx^2);          % (a.u.)^-2
%**************potential set up*****************************
V=zeros(N,1);
% Potential in Hr
for j=2:(N-1)
   V(j) = interp1(xscaleI,VI,xscale(j))/Hr;
end
V(1)=20;              %boundary condition
V(N)=20;              %boundary condition
%****************density of state     ****************
Do1 = g1*md1*m0/3.1415/(hb)^2/6.24146E11; %density of state (#/eV/cm2)
Do2 = g2*md2*m0/3.1415/(hb)^2/6.24146E11; %density of state (#/eV/cm2)
Do3 = g3*md3*m0/3.1415/(hb)^2/6.24146E11; %density of state (#/eV/cm2)
%************* Schrodinger Equation ********************
H = zeros(N,N);
%  light hole
   for j=2:(N-1)
      H(j,j) = V(j)+2*dd/m1;
   end
   H(N,N)=V(N)+2*dd/m1;
   H(1,1)=V(1)+2*dd/m1;
   for j=2:N
       H(j-1,j) = -dd/m1;
       H(j,j-1) = -dd/m1;
   end
   [Y,D]=eig(H); % Eigen vectors(Y) and Eigen values(D)
   [lambda1,key1] =sort(diag(D)); % Sort of eigen vector
   Y1 = Y(:,key1);
   E1=lambda1*Hr;
%  heavy hole
H = zeros(N,N);
   for j=2:(N-1)
      H(j,j) = V(j)+2*dd/m2;
   end
   H(N,N)=V(N)+2*dd/m2;
   H(1,1)=V(1)+2*dd/m2;
   for j=2:N
       H(j-1,j) = -dd/m2;
       H(j,j-1) = -dd/m2;
   end
   [Y,D]=eig(H); % Eigen vectors(Y) and Eigen values(D)
   [lambda2,key2] =sort(diag(D)); % Sort of eigen vector
   Y2 = Y(:,key2);
   E2=lambda2*Hr;
%  split band
H = zeros(N,N);
   for j=2:(N-1)
      H(j,j) = V(j)+2*dd/m3;
   end
   H(N,N)=V(N)+2*dd/m3;
   H(1,1)=V(1)+2*dd/m3;
   for j=2:N
       H(j-1,j) = -dd/m3;
       H(j,j-1) = -dd/m3;
   end
   [Y,D]=eig(H); % Eigen vectors(Y) and Eigen values(D)
   [lambda3,key3] =sort(diag(D)); % Sort of eigen vector
   Y3 = Y(:,key3);
   E3=lambda3*Hr+deltaE; % sprit is added
%**************** Calculating hole     densities ******************
%p:hole density /cm2
   for j=1:N
   p1(j)=Do1*k*T*log(1+exp((Ef-E1(j))/k/T)); %hole den in heavy
   p2(j)=Do2*k*T*log(1+exp((Ef-E2(j))/k/T)); %hole den in light
   p3(j)=Do3*k*T*log(1+exp((Ef-E3(j))/k/T)); %hole den in sprit
   end
%YY:hole density at each valley /cm3
   for j=1:N
    for jj=1:N
    YY1(j,jj)=(Y1(j,jj))^2*p1(jj)/dx0;  %hole den. (heavy band in #/cc)
    YY2(j,jj)=(Y2(j,jj))^2*p2(jj)/dx0;  %hole den. (light band in #/cc)
    YY3(j,jj)=(Y3(j,jj))^2*p3(jj)/dx0;  %hole den. (sprit off band in #/cc)
    end
   end
%R:hole density /cm3
   R=zeros(N,1);
   for j=1:N
    for jj=1:N
     R(j) = R(j)+(YY1(j,jj)+YY2(j,jj)+YY3(j,jj));
    end
   end
