function [IntGrand ] = GkIntPoly(s,P,Q,k)


IntGrand = besselh(0,1,k*sqrt(s.^2+P.*s+Q)); %greens function 


end