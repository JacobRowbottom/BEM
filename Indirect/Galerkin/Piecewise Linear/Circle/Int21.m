function [Int ] = Int21(s,t,si,tj,bi,bj,k,h,L)
if si==0 && bi==L
    si=L;
end
if tj==0 && bj==L
    tj=L;
end
E=double(eulergamma);
IntGrandpre=ones(size(t)).*(1+(2i/pi)*(E+log(k/2)));


Ker = (besselh(0,1,2*k*sin(abs(s-t)/2))-(2i/pi)*log(2*sin(abs(s-t)/2)));
Ker(abs(s-t) <1e-15)=IntGrandpre(abs(s-t)<1e-15);

Int=Ker.*((si-s+h)/h).*((t-tj+h)/h);
end