function [Int ] = PolyInt22(s,t,k,hi,hj,si,tj,bi,bj,L,sxv,syv,txv,tyv,CosEdgeAngleS,SinEdgeAngleS,CLS,CosEdgeAngleT,SinEdgeAngleT,CLT)
if si==0 && bi==L
    si=L;
end
if tj==0 && bj==L
    tj=L;
end
E=double(eulergamma);
% IntGrandpre=ones(size(t)).*(1+(2i/pi)*(E-log(2))) ;
IntGrandpre=ones(size(t)).*(1+(2i/pi)*(E+log(k/2)));

sx=sxv+(s-CLS)*CosEdgeAngleS; %cartesian coordinates in terms of arclength
sy=syv+(s-CLS)*SinEdgeAngleS;

tx=txv+(t-CLT)*CosEdgeAngleT; %cartesian coordinates in terms of arclength
ty=tyv+(t-CLT)*SinEdgeAngleT;

D=sqrt((sx-tx).^2+(sy-ty).^2);

Ker = besselh(0,1,k*D) - (2i/pi)*log(D); 
Ker(D<1e-15)=IntGrandpre(D<1e-15); 

Int=Ker.*((si-s+hi)/hi).*((tj-t+hj)/hj);

end