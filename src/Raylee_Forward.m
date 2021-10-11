function [vvp,U,vsv,evv] = ForWard(vsv,vpv,rhov,h,fks,Nn)

modn = 1;

countr = 0;
for f=fks

    [kk, vpk, vgk, ev] = ... 
        raylee_lysmer(Nn,vsv,vpv,rhov,f,h,modn,0,0,0,0);

    countr = countr + 1;
    vvp(1,countr) = vpk;
    U(1,countr) = vgk;
    evv(countr,:,:) = ev;
    
end
