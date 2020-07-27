x_post=[
   -0.3825
   -0.3408
   -0.3653
   -0.3285
   -0.3334
   -0.3016
   -0.3040
   -0.2697
   -0.2697
   -0.2256
   -0.2231
   -0.1422
   -0.1324
    0.0220
   -0.0245
    0.1274
    0.1274
    0.1765
    0.2059
    0.2230];
for i=1:10
    Vg_post(i)=x_post(2*i-1);
    Vgi_post(i)=x_post(2*i);
    delvgpost=Vg_post(i)-Vgi_post(i);
end;
for i=1:10
    for j=1:280
        if Vgg(j)<=Vgi_post(i) && Vgg(j+1)>=Vgi_post(i)
            psi_post(i)=(psi(j)+psi(j+1))/2;
            break;
        end
    end
end
            

