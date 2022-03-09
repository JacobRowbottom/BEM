function [IntGrand ] = GkIntPolyModD(si,k,xv,yv,CosEdgeAngle,SinEdgeAngle,CL,stepsize,xa,xb,ya,yb,a,b)

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
    
%     s1=0.5*(si(j)-a)*xg'+0.5*(si(j)+a);
%     w1=0.5*(si(j)-a)*wg;
%     
%     s2=0.5*(b-si(j))*xg'+0.5*(b+si(j));
%     w2=0.5*(b-si(j))*wg;

     st=0.5*(b-a)*xg'+0.5*(b+a);
     wt=0.5*(b-a)*wg;
    
    %IntGrand(j,1) = w1*GkIntPolyMod(s1,P(j),Q(j),k)+w2*GkIntPolyMod(s2,P(j),Q(j),k); %greens function 
    IntGrand(j,1) = wt*GkIntPolyMod(st,P(j),Q(j),k); %greens function 
end


end