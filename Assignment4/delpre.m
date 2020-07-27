x_pre=[    -0.3491
   -0.2195
   -0.3446
   -0.2106
   -0.3357
   -0.1882
   -0.3312
   -0.1748
   -0.2910
   -0.1525
   -0.2597
   -0.1123
   -0.1614
   -0.0542
   -0.0095
    0.0486
    0.1290
    0.1558
    0.2273
    0.2542
];
for i=1:10
    Vgi_pre(i)=x_pre(2*i-1);
    Vg_pre(i)=x_pre(2*i);
    delvgpre=Vg_pre(i)-Vgi_pre(i);
end;
for i=1:10
    for j=1:280
        if Vgg(j)<=Vgi_pre(i) && Vgg(j+1)>=Vgi_pre(i)
            psi_pre(i)=(psi(j)+psi(j+1))/2;
            break;
        end
    end
end