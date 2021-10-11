function [d] =  poisfun(a,b,c)

po = @(vp,vs) ((vp./vs)^2-2)./(2*(vp./vs)^2 - 2);
vs = @(vp,po) vp.*((0.5-po)/(1-po))^.5;
vp = @(vs,po) vs./((0.5-po)/(1-po))^0.5;

if c == 1;
    d = po(a,b);
elseif c ==2
    d = vs(a,b);
elseif c == 3
    d = vp(a,b);
end
    