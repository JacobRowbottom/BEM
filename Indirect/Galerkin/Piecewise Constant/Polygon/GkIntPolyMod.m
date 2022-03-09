function [IntGrand ] = GkIntPolyMod(s,P,Q,k)
E=double(eulergamma);

IntGrandpre=ones(size(s)).*(1+(2i/pi)*(E-log(2)));

IntGrand= besselh(0,1,k*sqrt(s.^2+P*s+Q))-(2*1i/pi)*log(sqrt(s.^2+P*s+Q));
 
IntGrand(sqrt(s.^2+P*s+Q)<1e-15)=IntGrandpre(sqrt(s.^2+P*s+Q)<1e-15); %New method 

% IntGrand = besselh(0,1,k*sqrt(s.^2+P*s+Q))-(2*1i/pi)*log(sqrt(s.^2+P*s+Q));  %Old method 

end