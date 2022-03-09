function [ IntGrand ] = HelmLinSolP(t,tj,xv,yv,CL,bj,L,CosEdgeAngle,SinEdgeAngle,k,h,x1,x2)
if tj==0 && bj==L
    tj=L;
end

xt=xv+(t-CL)*CosEdgeAngle;
yt=yv+(t-CL)*SinEdgeAngle;

D = sqrt((x1-xt).^2 + (x2 - yt).^2);

 IntGrand = besselh(0,1,k*D).*((tj-t+h)/h); %greens function 


end