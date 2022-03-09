function [IntGrand ] = GkPolyDM(t,tj,P,Q,k,stepsize) %ti is collocation, tj is t at point j
E=double(eulergamma);
IntGrandpre=ones(size(t)).*(1+(2i/pi)*(E-log(2)));

IntGrand = (besselh(0,1,k*sqrt(t.^2+P.*t+Q))-(2i/pi)*log(sqrt(t.^2+P.*t+Q))).*((tj-t+stepsize)/stepsize); %greens function 
 
IntGrand(sqrt(t.^2+P.*t+Q)<1e-15)=IntGrandpre(sqrt(t.^2+P.*t+Q)<1e-15);
end