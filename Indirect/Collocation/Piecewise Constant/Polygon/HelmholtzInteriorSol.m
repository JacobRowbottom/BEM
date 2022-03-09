function [ IntGrand ] = HelmholtzInteriorSol( s,xv,yv,CL,CosEdgeAngle,SinEdgeAngle,k,x1,x2)


xs=xv+(s-CL)*CosEdgeAngle;
ys=yv+(s-CL)*SinEdgeAngle;

IntGrand = (-1i/4)*besselh(0,1,k*sqrt((x1-xs).^2 + (x2 - ys).^2)); %greens function 


end
