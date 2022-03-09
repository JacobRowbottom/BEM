function [IntGrand ] = G0tNum22(s,t,tj,bj,h,L) %ti is collocation, tj is t at point j
if tj==0 && bj==L
    tj=L;
end
IntGrandpre=zeros(size(t));

Ker= log(2.*sin(abs(s-t)/2));
Ker(abs(s-t) <1e-15)=IntGrandpre(abs(s-t)<1e-15); %New method 

IntGrand=Ker.*((t-s)/h).*((tj-t+h)/h);
end