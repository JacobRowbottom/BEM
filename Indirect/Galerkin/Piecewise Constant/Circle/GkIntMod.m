function [IntGrand ] = GkIntMod(s,si,k)


IntGrand = besselh(0,1,2*k*sin(abs(si-s)/2))-(2*1i/pi)*log(2*sin(abs(si-s)/2)); %greens function 


end