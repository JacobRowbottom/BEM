function [IntGrand ] = GkIntPolyM(t,tj,P,Q,k,stepsize)
E=double(eulergamma);
IntGrandpre=ones(size(t)).*(1+(2i/pi)*(E-log(2)));

IntGrand = besselh(0,1,k*sqrt(t.^2+P.*t+Q)).*((tj-t+stepsize)/stepsize); %greens function 
 
IntGrand(sqrt(t.^2+P.*t+Q)<1e-15)=IntGrandpre(sqrt(t.^2+P.*t+Q)<1e-15);

end
