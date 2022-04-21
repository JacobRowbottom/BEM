clear 

SL=zeros(1,4);
nEdge=zeros(1,4);
h=zeros(1,4);
CosEdgeAngle=zeros(1,4);
SinEdgeAngle=zeros(1,4);
CL=zeros(1,4);
cEdge=zeros(1,4);

 xv(1)=0; yv(1)=0; %coordinates of vertices 
 xv(2)=1; yv(2)=0;
 xv(3)=1; yv(3)=1;
 xv(4)=0; yv(4)=1;

% xv = [0 1 1 0.5 0.5 0]; % L- shape
% yv = [0 0 1 1 0.5 0.5];

nV=length(xv); %no. of vertices of polygons 

NA=64; %number of elements
k=1;   %wavenumber  

for j=1:nV
    jp = j+1;
    if (jp > nV)
        jp=jp-nV;
    end
    SL(j)=sqrt((xv(j)-xv(jp))^2+(yv(j)-yv(jp))^2); %side lengths 
    CL(j)=sum(SL(1:j-1)); %cummalative length 
    CosEdgeAngle(j)=(xv(jp)-xv(j))/SL(j);
    SinEdgeAngle(j)=(yv(jp)-yv(j))/SL(j);
    
end

L=sum(SL); %perimeter of polygon
Ds = L/NA; %average length of element/ 

for j=1:nV
    nEdge(j) = round(SL(j)/Ds); %Number of elements on each side
    cEdge(j)=sum(nEdge(1:j-1));   % cummaltive number of elements on edge
    h(j)=SL(j)/nEdge(j); %stepsize on each side
    
    for je=1:nEdge(j)
        
        a(je+cEdge(j))=CL(j)+(je-1)*h(j); %element end points 
        b(je+cEdge(j))=CL(j)+je*h(j);
        
        xa(je+cEdge(j))=xv(j)+h(j)*CosEdgeAngle(j)*(je-1);
        ya(je+cEdge(j))=yv(j)+h(j)*SinEdgeAngle(j)*(je-1); %coordinates of elements

        xb(je+cEdge(j))=xv(j)+h(j)*CosEdgeAngle(j)*(je);
        yb(je+cEdge(j))=yv(j)+h(j)*SinEdgeAngle(j)*(je);

    end
end

s=a; xs=xa; ys=ya;
tj=a;

N=sum(nEdge);
for i=1:nV
     if i==1
       f(1+cEdge(i):cEdge(i)+nEdge(i),1) = cos(k*(xs(1+cEdge(i):cEdge(i)+nEdge(i))))'; % boundary conditions for square 
    elseif i==2
       f(1+cEdge(i):cEdge(i)+nEdge(i),1)=cos(k)*ones(nEdge(i),1); 
    elseif i==3
        f(1+cEdge(i):cEdge(i)+nEdge(i),1) = cos(k*(xs(1+cEdge(i):cEdge(i)+nEdge(i))))';
    else
        f(1+cEdge(i):cEdge(i)+nEdge(i),1)=ones(nEdge(i),1); 
     end
%         if i==6
%             f(1+cEdge(i):cEdge(i)+nEdge(i),1)=cos(k*(xs(1+cEdge(i):cEdge(i)+nEdge(i))))'; % vector of L shape
%         else 
%             f(1+cEdge(i):cEdge(i)+nEdge(i),1)=0;
%         end 
end

for i=1:nV
    for iv=1:nEdge(i)
         it=iv+cEdge(i);
         for j=1:nV
            for jv=1:nEdge(j)
                jt=jv+cEdge(j);
                jtp=jt+1;
                if jtp>NA
                    jtp=jtp-NA;
                end
                if it==jt
                    MM(it,jt)=(-1i/4)*(integral(@(t)PolyIntDM(t,k,h(j),xs(jt),ys(jt),tj(jtp),b(jt),L,xv(j),yv(j),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(jt),b(jt)) + (2i/pi)*(integral(@(t)G0tInt(s(it),t,tj(jtp),b(jt),h(j),L),a(jt),b(jt)) + G0ExactPolyM(a(jt),b(jt),s(it),tj(jtp),b(jt),h(j),L)));
                    MP(it,jt)=(-1i/4)*(integral(@(t)PolyIntDP(t,k,h(j),xs(jt),ys(jt),tj(jt),b(jt),L,xv(j),yv(j),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(jt),b(jt)) + (2i/pi)*(integral(@(t)G0tInt(s(it),t,tj(jt),b(jt),h(j),L),a(jt),b(jt)) + G0ExactPolyM(a(jt),b(jt),s(it),tj(jt),b(jt),h(j),L)));
                else
                    MM(it,jt) = (-1i/4)*integral(@(t)PolyIntM(t,k,h(j),xs(it),ys(it),tj(jtp),b(jt),L,xv(j),yv(j),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(jt),b(jt));
                    MP(it,jt) = (-1i/4)*integral(@(t)PolyIntP(t,k,h(j),xs(it),ys(it),tj(jt),b(jt),L,xv(j),yv(j),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(jt),b(jt));
                end
            end
         end

    end
end

MFinal = MM+MP(1:end,[2:end,1]);
sigma =MFinal\f;            
 
sol=1; %pick from sol=0 to solve at 1 point or sol=1 to solve on many points of a grid
for i=1:nV
    for iv=1:nEdge(i)
        it=iv+cEdge(i); 
        EdgeIndex(it,1)=i;
    end
end

if sol==0
    
x1 = 0.75; % x cartesian coordinates
x2 = 0.25;
uInter=zeros(1,N);
   for i =1:nV          
    for iv=1:nEdge(i)
            it=iv+cEdge(i);
            itp=it+1;
            if itp>NA
               itp=itp-NA;
            end
            uInter(it) = (-1i/4)*(integral(@(t)HelmLinSolM(t,tj(itp),xv(i),yv(i),CL(i),b(it),L,CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),x1,x2),a(it),b(it))+integral(@(t)HelmLinSolP(t,tj(itp),xv(EdgeIndex(itp)),yv(EdgeIndex(itp)),CL(EdgeIndex(itp)),b(itp),L,CosEdgeAngle(EdgeIndex(itp)),SinEdgeAngle(EdgeIndex(itp)),k,h(EdgeIndex(itp)),x1,x2),a(itp),b(itp)));
            
    end
   end
   u = uInter*sigma

   ue = cos(k*x1); % Exact solution 

elseif sol==1
    
    nplot=25;
    v1=linspace(0,1,nplot); %gird for square
    v2=linspace(0,1,nplot);
%     
%     v1=linspace(-1,1,nplot); %grid for L shape
%     v2=linspace(-1,1,nplot);

    [V1,V2] =meshgrid(v1,v2);
    
    for j=1:nplot;
        j
         for m=1:nplot;
             for i =1:nV  
                 for iv=1:nEdge(i)
                     it=iv+cEdge(i);
                     itp=it+1;
                     if itp>NA
                        itp=itp-NA;
                     end
%                      uInter(it) = (-1i/4)*(integral(@(t) HelmLinSolP(t,s(it),xv(im),yv(im),CL(im),CosEdgeAngle(im),SinEdgeAngle(im),k,h(im),v1(j),v2(m),b(itm),L),a(itm),b(itm)) + integral(@(t) HelmLinSolM(t,s(it),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),v1(j),v2(m)),a(it),b(it)));
                      uInter(it) = (-1i/4)*(integral(@(t)HelmLinSolM(t,tj(itp),xv(i),yv(i),CL(i),b(it),L,CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),v1(j),v2(m)),a(it),b(it))+integral(@(t)HelmLinSolP(t,tj(itp),xv(EdgeIndex(itp)),yv(EdgeIndex(itp)),CL(EdgeIndex(itp)),b(itp),L,CosEdgeAngle(EdgeIndex(itp)),SinEdgeAngle(EdgeIndex(itp)),k,h(EdgeIndex(itp)),v1(j),v2(m)),a(itp),b(itp)));


                 end
             end
             u(j,m)= uInter*sigma;
             ue(j,m) = cos(k*v1(j));
          end
    end
     figure
     mesh(V2,V1,real(u))
     figure
     mesh(V2,V1,real(ue))
end                 
                           
                 
