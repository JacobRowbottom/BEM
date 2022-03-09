function [ val ] = G0Int(a,b,s) %%%% Exact answer for the integral of  ln(abs(s-t)) for circle case. 
if s==b
valb= -1i.*polylog(2,1);
vala=((s - a).*log(1-exp(1i.*(s - a))) - (s-a).*log(-2*sin(0.5.*(s-a))) - 1i*(0.25.*(s-a).^2 + polylog(2,exp(1i*(s-a))))).*(0.5-0.5*sign(s-a)) + (((s -a).*log(1-exp(1i*(s-a))) - (s-a).*log(2*sin(0.5*(s -a))) - 1i*(0.25*(s -a).^2 + polylog(2,exp(1i*(s -a))))).*(1 + sign(s -a)))/2;

elseif s==a 
valb=((s - b).*log(1-exp(1i.*(s - b))) - (s-b).*log(-2*sin(0.5.*(s-b))) - 1i*(0.25.*(s-b).^2 + polylog(2,exp(1i*(s-b))))).*(0.5-0.5*sign(s-b)) + (((s -b).*log(1-exp(1i*(s-b))) - (s-b).*log(2*sin(0.5*(s -b))) - 1i*(0.25*(s -b).^2 + polylog(2,exp(1i*(s -b))))).*(1 + sign(s -b)))/2;
vala= -1i.*polylog(2,1);
   
else 
valb=((s - b).*log(1-exp(1i.*(s - b))) - (s-b).*log(-2*sin(0.5.*(s-b))) - 1i*(0.25.*(s-b).^2 + polylog(2,exp(1i*(s-b))))).*(0.5-0.5*sign(s-b)) + (((s -b).*log(1-exp(1i*(s-b))) - (s-b).*log(2*sin(0.5*(s -b))) - 1i*(0.25*(s -b).^2 + polylog(2,exp(1i*(s -b))))).*(1 + sign(s -b)))/2;
vala=((s - a).*log(1-exp(1i.*(s - a))) - (s-a).*log(-2*sin(0.5.*(s-a))) - 1i*(0.25.*(s-a).^2 + polylog(2,exp(1i*(s-a))))).*(0.5-0.5*sign(s-a)) + (((s -a).*log(1-exp(1i*(s-a))) - (s-a).*log(2*sin(0.5*(s -a))) - 1i*(0.25*(s -a).^2 + polylog(2,exp(1i*(s -a))))).*(1 + sign(s -a)))/2;
    
end
val=valb-vala;

end 
