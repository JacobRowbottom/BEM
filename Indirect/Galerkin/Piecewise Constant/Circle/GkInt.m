function [IntGrand ] = GkInt(s,si,k)


IntGrand = besselh(0,1,2*k*sin(abs(si-s)/2)); %greens function 


end