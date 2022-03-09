function [IntGrand ] = GkIntPolyGal(si,k,xv,yv,CosEdgeAngle,SinEdgeAngle,CL,stepsize,xa,xb,ya,yb,a,b)

xi=xv+(si-CL)*CosEdgeAngle;
yi=yv+(si-CL)*SinEdgeAngle;
 
da=sqrt((xi-xa).^2+(yi-ya).^2);
db=sqrt((xi-xb).^2+(yi-yb).^2);
P=((db.^2-da.^2 -stepsize.^2)/stepsize)-2*a;
Q=da.^2-((db.^2-da.^2 -stepsize.^2)*a/stepsize)+a.^2;

N=length(si);
nGaussJ=17;
[xg, wg]=grule(nGaussJ);

for j=1:N
    
    sj=0.5*(b-a)*xg'+0.5*(b+a);
    wj=0.5*(b-a)*wg;
    
    IntGrand(j,1) = wj*GkIntPoly(sj,P(j),Q(j),k);%greens function 
end

end