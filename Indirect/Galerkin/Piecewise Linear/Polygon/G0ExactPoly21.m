function [ Int ] = G0ExactPoly21(a,b,t,si,tj,bj,hi,hj,L)%s is the arc length 
if tj==0 && bj==L
    tj=L;
end
if si==0 && b==L
    si=L;
end

if t==b 
valb=0;
vala=(a-t).*(log(abs(a-t))-1);
elseif t==a
valb=(b-t).*(log(abs(b-t))-1);
vala=0;
else 
valb=(b-t).*(log(abs(b-t))-1);
vala=(a-t).*(log(abs(a-t))-1);  
end

val=valb-vala;

Int=((si-t+hi)/hi).*val.*((t-tj+hj)/hj);

end