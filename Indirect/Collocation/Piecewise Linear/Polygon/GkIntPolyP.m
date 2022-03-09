function [IntGrand ] = GkIntPolyP(t,tj,b,P,Q,k,stepsize,L)
E=double(eulergamma);
if tj==0 && b==L
    tj=L;
end

IntGrandpre=ones(size(t)).*(1+(2i/pi)*(E-log(2)));

IntGrand = besselh(0,1,k*sqrt(t.^2+P.*t+Q)).*((t-tj+stepsize)/stepsize); %greens function 

IntGrand(sqrt(t.^2+P.*t+Q)<1e-15)=IntGrandpre(sqrt(t.^2+P.*t+Q)<1e-15);

end