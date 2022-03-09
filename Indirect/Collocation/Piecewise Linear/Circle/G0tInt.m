function [IntGrand ] = G0tInt(t,ti,tj,h) %ti is collocation, tj is t at point j
% if tj==0 && bj==L
%     tj=L;
% end
IntGrandpre=zeros(size(t));

IntGrand= log(2*sin(abs(ti-t)/2)).*((t-tj)/h);
 
IntGrand(abs(ti-t)/2 <1e-15)=IntGrandpre(abs(ti-t)/2<1e-15); %New method 
end