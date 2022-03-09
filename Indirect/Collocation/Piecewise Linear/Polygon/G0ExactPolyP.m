function [ Int ] = G0ExactPolyP(a,b,si,tj,bj,h,L)%s is the arc length 
if tj==0 && bj==L
    tj=L;
end

if si==b 
valb=0;
vala=(a-si).*(log(abs(a-si))-1);
elseif si==a
valb=(b-si).*(log(abs(b-si))-1);
vala=0;
else 
valb=(b-si).*(log(abs(b-si))-1);
vala=(a-si).*(log(abs(a-si))-1);  
end

val=valb-vala;

Int=((tj-si+h)/h).*val;

end