function [IntGrand ] = GkIntLin(t,ti,tj,k,h) %ti is collocation, tj is t at point j

IntGrand = besselh(0,1,2*k*sin(abs(ti-t)/2)).*((tj+h-t)/h); %greens function 

end