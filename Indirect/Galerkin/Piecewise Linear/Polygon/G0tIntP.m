function [IntGrand ] = G0tIntP(t,si,tj,bj,hi,hj,L) %ti is collocation, tj is t at point j
if tj==0 && bj==L
    tj=L;
end
IntGrandpre=zeros(size(t));

Ker= log(abs(si-t));
Ker(abs(si-t) <1e-15)=IntGrandpre(abs(si-t)<1e-15); %New method 

IntGrand=((si-t)/hi).*Ker.*((t-tj+hj)/hj);
end