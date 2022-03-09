function [ Int ] = HelmLinSolP( t,tj,xv,yv,CL,CosEdgeAngle,SinEdgeAngle,k,h,x1,x2,b,L)

if tj==0 && b==L
    tj=L;
end

xt=xv+(t-CL)*CosEdgeAngle;
yt=yv+(t-CL)*SinEdgeAngle;

D=sqrt((x1-xt).^2 + (x2-yt).^2);

Int = besselh(0,1,k*D).*((t-tj+h)/h); %greens function 


end