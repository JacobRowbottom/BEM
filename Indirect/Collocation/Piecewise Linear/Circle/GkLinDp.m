function [IntGrand ] = GkLinDp(t,ti,tj,k,h) %ti is collocation, tj is t at point j
E=double(eulergamma);
IntGrandpre=ones(size(t)).*(1+(2i/pi)*(E-log(2)));

IntGrand= (besselh(0,1,2*k*sin(abs(ti-t)/2))-(2i/pi)*log(2*sin(abs(ti-t)/2))).*((t-tj+h)/h);
 
IntGrand(abs(ti-t)/2 <1e-15)=IntGrandpre(abs(ti-t)/2<1e-15); %New method
end