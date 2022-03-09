function [Int ] = PolyIntDP(t,k,h,sx,sy,tj,bj,L,txv,tyv,CosEdgeAngleT,SinEdgeAngleT,CLT)
if tj==0 && bj==L
    tj=L;
end
E=double(eulergamma);
IntGrandpre=ones(size(t)).*(1+(2i/pi)*(E+log(k/2)));

tx=txv+(t-CLT)*CosEdgeAngleT; %cartesian coordinates in terms of arclength
ty=tyv+(t-CLT)*SinEdgeAngleT;

D=sqrt((sx-tx).^2+(sy-ty).^2);

Ker = besselh(0,1,k*D) - (2i/pi)*log(D); 
Ker(D<1e-15)=IntGrandpre(D<1e-15); 

Int=Ker.*((tj-t+h)/h);

end