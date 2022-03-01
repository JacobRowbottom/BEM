function [Ker] = dGkdn(k,xi,yi,nxq,nyq,x0,y0)
%x0, y0 point source coords 
% xi, yi coll pt coords 

D=sqrt((xi-x0).^2+(yi-y0).^2);
Bessel = (-k).*besselh(1,1,k.*D);

Ker = (1./D).*(nxq.*Bessel.*(xi-x0) +  nyq.*Bessel.*(yi-y0));
end