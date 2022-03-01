NP=25;

uInter=zeros(NP,NP);
Gkfun=zeros(NP,NP);
VSLPint = zeros(M,1);  
KDLPint = zeros(M,1);

% IPx=0.5;
% IPy=0.5;

IPx=linspace(0,1,NP+2); %gird for square
IPy=linspace(0,1,NP+2);

IPx = IPx(2:end-1);
IPy = IPy(2:end-1);

[IPX,IPY] =meshgrid(IPx,IPy);
for ipx=1:NP
        ipx
         for ipy=1:NP
             
            for i=1:NVert
                for iv=1:nEdge(i)
                    it=iv+cEdge(i);
                    
                    VSLPint(it)=(-1i/4)*integral(@(q)VSLPIntPoly(q,k,IPx(ipx),IPy(ipy),xv(i),yv(i),CosEdgeAngle(i),SinEdgeAngle(i),CL(i)),a(it),b(it));     
                    KDLPint(it)=(-1i/4)*integral(@(q)KIntPoly(q,k,IPx(ipx),IPy(ipy),xv(i),yv(i),CosEdgeAngle(i),SinEdgeAngle(i),CL(i),nxq(it),nyq(it)),a(it),b(it));                   
                end
            end                                 
             uInter(ipx,ipy) = ((vBC.')*VSLPint-(vBound.')*KDLPint);

             if IPx(ipx)<0.5 && IPy(ipy)>0.5 
                 uInter(ipx,ipy)=0;
                 Gkfun(ipx,ipy)=0;
             else
                 Gkfun(ipx,ipy) = -(1i/4)*Gk(k,IPx(ipx),IPy(ipy),PSx,PSy);
             end
         end

end
uPS = uInter + Gkfun;
figure
mesh(IPx,IPX,real(uPS))
