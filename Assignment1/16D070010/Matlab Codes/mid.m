Vg=-1.5;
epsi=12*8.89*10^-14;
KbTq=0.025;
ni=10^10;
q=1.6*10^-19;
Cox=4*8.89*10^-14/(2*10^-7);
co=sqrt(2*epsi*KbTq*q*ni)/Cox;
syms psis;
ans=vpasolve(-Vg+psis-co*((exp(psis/KbTq)+exp(-psis/KbTq)-2)^0.5)==0, psis);
ans1=(co*Cox)*((exp(ans/KbTq)+exp(-ans/KbTq)-2)^0.5);


