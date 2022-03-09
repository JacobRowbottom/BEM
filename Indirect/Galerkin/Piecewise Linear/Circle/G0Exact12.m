function [ Int ] = G0Exact12(a,b,t,si,bi,tj,bj,h,L) %%%% Exact answer for the integral of  ln(abs(s-t)) for circle case. 
if tj==0 && bj==L
    tj=L;
end
if si==0 && bi==L
    si=L;
end
if t==b
valb= -1i.*polylog(2,1);
vala=((t - a).*log(1-exp(1i.*(t - a))) - (t-a).*log(-2*sin(0.5.*(t-a))) - 1i*(0.25.*(t-a).^2 + polylog(2,exp(1i*(t-a))))).*(0.5-0.5*sign(t-a)) + (((t -a).*log(1-exp(1i*(t-a))) - (t-a).*log(2*sin(0.5*(t -a))) - 1i*(0.25*(t -a).^2 + polylog(2,exp(1i*(t -a))))).*(1 + sign(t -a)))/2;

elseif t==a
valb=((t - b).*log(1-exp(1i.*(t - b))) - (t-b).*log(-2*sin(0.5.*(t-b))) - 1i*(0.25.*(t-b).^2 + polylog(2,exp(1i*(t-b))))).*(0.5-0.5*sign(t-b)) + (((t -b).*log(1-exp(1i*(t-b))) - (t-b).*log(2*sin(0.5*(t -b))) - 1i*(0.25*(t -b).^2 + polylog(2,exp(1i*(t -b))))).*(1 + sign(t -b)))/2;
vala= -1i.*polylog(2,1);
   
else 
valb=((t - b).*log(1-exp(1i.*(t - b))) - (t-b).*log(-2*sin(0.5.*(t-b))) - 1i*(0.25.*(t-b).^2 + polylog(2,exp(1i*(t-b))))).*(0.5-0.5*sign(t-b)) + (((t -b).*log(1-exp(1i*(t-b))) - (t-b).*log(2*sin(0.5*(t -b))) - 1i*(0.25*(t -b).^2 + polylog(2,exp(1i*(t -b))))).*(1 + sign(t -b)))/2;
vala=((t - a).*log(1-exp(1i.*(t - a))) - (t-a).*log(-2*sin(0.5.*(t-a))) - 1i*(0.25.*(t-a).^2 + polylog(2,exp(1i*(t-a))))).*(0.5-0.5*sign(t-a)) + (((t -a).*log(1-exp(1i*(t-a))) - (t-a).*log(2*sin(0.5*(t -a))) - 1i*(0.25*(t -a).^2 + polylog(2,exp(1i*(t -a))))).*(1 + sign(t -a)))/2;
    
end
val=valb-vala;

Int=((h-si+t)/h).*val.*((tj-t+h)/h);


end 