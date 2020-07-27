function [xscale,E1,E2,Y1,Y2,YY1,YY2,R]=shelec(xscaleI,VI,xstart,xend,N,Ef,T);
%***********************************************************
% 1-dimensions Shrodinger equations
% for (100) electorn ver 1.0.0
% "input"
% xscaleI:cordinate of volatge (cm)
% VI: voltage (V)
% (xstrat,xend): analysis region (cm)
% N: analysis mesh number
% Ef: Fermi level
% T : temperature
% "output"
% xscale : new scale (cm)
% E1-2 (i):energy level of light, heavy on i-band
% YY1-2   (j,i):density of light, heavy on i-band at j-point
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
Nv = 1.02E19;      % effective density of states of valence band
Nc = 2.8E19;       % effective density of states of conduction band 
m0 = 9.109534E-28; % electron mass in g
m1 = 0.916;        % Electron mass in z direction for lower(3,6)
m2 = 0.19;         % Electron mass in z direction for higher(1,4,2,5)
md1 = 0.19;        % Density of state mass low energy valley (3,6)
md2 = 0.417;       % Density of state mass high energy valley (1,4,2,5)
g1 = 2;            % degeneracy for low energy valley (3,6)
g2 = 4;            % degeneracy for high energy valley (1,4,2,5) 
%****************      set up         ****************
Ef=-Ef;
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
%****************CALCULATED PARAMETERS****************
Do1 = g1*md1*m0/3.1415/(hb)^2/6.24146E11; %density of state (#/eV/cm2)
Do2 = g2*md2*m0/3.1415/(hb)^2/6.24146E11; %density of state (#/eV/cm2)
%************* Schrodinger Equation ********************
H = zeros(N,N);
%  light electorn
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
%  heavy electorn
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
%**************** Calculating electron densities ******************
%p:electron density /cm2
   for j=1:N
   ne1(j)=Do1*k*T*log(1+exp((Ef-E1(j))/k/T));
   ne2(j)=Do2*k*T*log(1+exp((Ef-E2(j))/k/T));
   end
%YY:electron density at each valley /cm3
   for j=1:N
    for jj=1:N
    YY1(j,jj)=(Y1(j,jj))^2*ne1(jj)/dx0;  %elec. den. (low ene. valley in #/cc)
    YY2(j,jj)=(Y2(j,jj))^2*ne2(jj)/dx0;  %elec. den. (high ene. valley in #/cc)
    end
   end
%R:electron density /cm3
   R = zeros (N,1);
   for j=1:N
    for jj=1:N 
      R(j) =R(j)+(YY1(j,jj)+YY2(j,jj));
    end
   end
