function [Ker] = Gk(k,xi,yi,x0,y0)
%x0, y0 point source coords 
% xi, yi coll pt coords 

D=sqrt((xi-x0).^2+(yi-y0).^2);
Ker = besselh(0,1,k.*D);


end