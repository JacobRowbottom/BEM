function [ val ] = G0IntPoly(a,b,s,P,Q,L)

%valb= (0.5*sqrt((4*Q - A^2))*atan2((A + 2*b),(sqrt((4*B - A^2)))) + 0.25*(A+2*b)*log(b*(A+b)+ B) - b);
%vala= (0.5*sqrt((4*B - A^2))*atan2((A + 2*a),(sqrt((4*B - A^2)))) + 0.25*(A+2*a)*log(a*(A+a)+ B) - a);

%same as above.
if s==0 && b==L
    s=L;
end

if s==b 
valb=0;
vala= a*log((a^2 + P*a + Q)^(1/2)) - a + (P*log(a^2 + P*a + Q))/4 + 2*atan2((P*(- P^2/16 + Q/4)^(1/2))+ (2*a*(- P^2/16 + Q/4)^(1/2)),(- P^2/4 + Q))*(- P^2/16 + Q/4)^(1/2);
elseif s==a
valb= b*log((b^2 + P*b + Q)^(1/2)) - b + (P*log(b^2 + P*b + Q))/4 + 2*atan2((P*(- P^2/16 + Q/4)^(1/2))+ (2*b*(- P^2/16 + Q/4)^(1/2)),(- P^2/4 + Q))*(- P^2/16 + Q/4)^(1/2);
vala=0;
else 
valb= b*log((b^2 + P*b + Q)^(1/2)) - b + (P*log(b^2 + P*b + Q))/4 + 2*atan2((P*(- P^2/16 + Q/4)^(1/2))+ (2*b*(- P^2/16 + Q/4)^(1/2)),(- P^2/4 + Q))*(- P^2/16 + Q/4)^(1/2);
vala= a*log((a^2 + P*a + Q)^(1/2)) - a + (P*log(a^2 + P*a + Q))/4 + 2*atan2((P*(- P^2/16 + Q/4)^(1/2))+ (2*a*(- P^2/16 + Q/4)^(1/2)),(- P^2/4 + Q))*(- P^2/16 + Q/4)^(1/2);
  
end

val=valb-vala;

% valb= b*log((b^2 + P*b + Q)^(1/2)) - b + (P*log(b^2 + P*b + Q))/4 + 2*atan2((P*(- P^2/16 + Q/4)^(1/2))+ (2*b*(- P^2/16 + Q/4)^(1/2)),(- P^2/4 + Q))*(- P^2/16 + Q/4)^(1/2);
% vala= a*log((a^2 + P*a + Q)^(1/2)) - a + (P*log(a^2 + P*a + Q))/4 + 2*atan2((P*(- P^2/16 + Q/4)^(1/2))+ (2*a*(- P^2/16 + Q/4)^(1/2)),(- P^2/4 + Q))*(- P^2/16 + Q/4)^(1/2);
%  
% 
% val=valb-vala;
 
end