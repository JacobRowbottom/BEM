function [IntGrand ] = G0tIntM(s,t,tj,bj,hi,hj,L) %ti is collocation, tj is t at point j
if tj==0 && bj==L
    tj=L;
end
IntGrandpre=zeros(size(t));

Ker= log(abs(s-t));
Ker(abs(s-t) <1e-15)=IntGrandpre(abs(s-t)<1e-15); %New method 

IntGrand=((s-t)/hi).*Ker.*((tj-t+hj)/hj);
end