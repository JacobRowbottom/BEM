function [ val ] = G0IntPolyD(a,b,s)

% valb= (-2*b + b*log(s - b) - 2*s*log(b-s) + b*log(b-s) + b*(log(s - b) - log(b-s))*sign(s - b))/2;
% vala= (-2*a + a*log(s - a) - 2*s*log(a-s) + a*log(a-s) + a*(log(s - a) - log(a-s))*sign(s - a))/2;
%  
%valb= (0.5)*(b*sign(s-b)*log(s-b) +b*sign(b-s)*log(b-s)+b*log(s-b)+b*log(b-s)-2*s*log(b-s)-2*b);
%vala= (0.5)*(a*sign(s-a)*log(s-a) +a*sign(a-s)*log(a-s)+a*log(s-a)+a*log(a-s)-2*s*log(a-s)-2*a);
if s==b
valb=0; %valb = -1/s      potentially ???
vala=(a-s).*(log(abs(a-s))-1);
%vala= -stepsize.*(log(abs(stepsize))-1);
elseif s==a
valb=(b-s).*(log(abs(b-s))-1);
vala=0;
else 
valb=(b-s).*(log(abs(b-s))-1);
vala=(a-s).*(log(abs(a-s))-1);  
end

val=valb-vala;
 
end

