NP=1;
IPx=0.5;
IPy=0.5;

VSLPint = zeros(M,1);  
KDLPint = zeros(M,1);
uInter=zeros(NP,NP);
uExact=zeros(NP,NP);

% IPx=linspace(0,1,NP+2); %gird for square
% IPy=linspace(0,1,NP+2);
% 
% IPx = IPx(2:end-1);
% IPy = IPy(2:end-1);

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
             uInter(ipx,ipy) = ((f.')*VSLPint-(u.')*KDLPint);
%              uExact(ipx,ipy) = (1/k)*(cot(k)*cos(k*IPx(ipx)) + sin(k*IPx(ipx))); %Unsure if this is correct solution for interior point

         end

end
% figure
% mesh(IPY,IPX,real(uInter))
% figure
% mesh(IPY,IPX,real(uExact))
