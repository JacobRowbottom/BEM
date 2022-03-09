function [ val ] = G0IntPolyD(a,b,s,L)%s is the arc length 
%function [ val ] = G0IntPolyD(a,b,s)
% valb= (-2*b + b*log(s - b) - 2*s*log(b-s) + b*log(b-s) + b*(log(s - b) - log(b-s))*sign(s - b))/2;
% vala= (-2*a + a*log(s - a) - 2*s*log(a-s) + a*log(a-s) + a*(log(s - a) - log(a-s))*sign(s - a))/2;


%%%%% integration of log(abs(t-ti)), ti is collocation point. 
if s==0 && b==L
    s=L;
end

if s==b 
valb=0;
vala=(a-s).*(log(abs(a-s))-1);
elseif s==a
valb=(b-s).*(log(abs(b-s))-1);
vala=0;
else 
valb=(b-s)*(log(abs(b-s))-1);
vala=(a-s)*(log(abs(a-s))-1);  
end

val=valb-vala;
 
end

