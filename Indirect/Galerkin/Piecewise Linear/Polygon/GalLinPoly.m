clear 
tic

SL=zeros(1,4);
nEdge=zeros(1,4);
h=zeros(1,4);
CosEdgeAngle=zeros(1,4);
SinEdgeAngle=zeros(1,4);
CL=zeros(1,4);
cEdge=zeros(1,4);
%coordinates of vertices
%Square 
 xv(1)=0; yv(1)=0;  
 xv(2)=1; yv(2)=0;
 xv(3)=1; yv(3)=1;
 xv(4)=0; yv(4)=1;
%Lshape
%  xv = [0 1 1 0.5 0.5 0];
%  yv = [0 0 1 1 0.5 0.5];

nV=length(xv); %no. of vertices of polygons 

NA=8; %number of elements
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
u=a;

N=sum(nEdge); %Total number of elements

for i=1:nV
    for iv=1:nEdge(i)
        it=iv+cEdge(i);
        itp=it+1;
        if itp>N
            itp=itp-N;
        end
     %Dirichlet boundary condition for each edge.
        if i==1
            fp(it,1)= integral(@(x) FCosPt(x,u(itp),b(it),k,h(i),L,CL(i),xv(i),CosEdgeAngle(i)),a(it),b(it));
            fm(it,1)= integral(@(x) FCosMt(x,u(it),b(it),k,h(i),L,CL(i),xv(i),CosEdgeAngle(i)),a(it),b(it));

        elseif i==2
            fp(it,1)= (cos(k)/h(i)).*((b(it).^2 - a(it).^2)/2 + (b(it)-a(it)).*(h(i)-u(itp)));
            fm(it,1)= (cos(k)/h(i)).*((u(it)+h(i)).*(b(it)-a(it)) + (a(it).^2 - b(it).^2)/2);
        elseif i==3
            fp(it,1)= integral(@(x) FCosPt(x,u(itp),b(it),k,h(i),L,CL(i),xv(i),CosEdgeAngle(i)),a(it),b(it));
            fm(it,1)= integral(@(x) FCosMt(x,u(it),b(it),k,h(i),L,CL(i),xv(i),CosEdgeAngle(i)),a(it),b(it));
        else
            if u(itp)==0 && b(it)==L
                u(itp)=L;
            end
            fp(it,1)= (1/h(i)).*((b(it).^2 - a(it).^2)/2 + (b(it)-a(it)).*(h(i)-u(itp))); 
            fm(it,1)= (1/h(i)).*((u(it)+h(i)).*(b(it)-a(it)) - ((b(it).^2 -a(it).^2)/2));
        end
%         if i==3     % vector of L shape
%             fp(it,1)= integral(@(x) FCosPt(x,u(itp),b(it),k,h(i),L,CL(i),xv(i),CosEdgeAngle(i)),a(it),b(it));
%             fm(it,1)= integral(@(x) FCosMt(x,u(it),b(it),k,h(i),L,CL(i),xv(i),CosEdgeAngle(i)),a(it),b(it)); 
%         else 
%             fp(it,1)=0;
%             fm(it,1)=0;
%         end
        f= fp([end 1: end-1]) + fm(1:end); %%% Calculate f
    end
end
u=a; 
%Solved as in Appendix B
MatrixGalLinPoly;

%1 indicates the basis function (t - tj-1)/h 
%2 indicates the basis function (tj+1 - t)/h 
%in Eqn (2.36) in thesis.
for i=1:nV
    im=i-1;
    if im<1
        im=im+nV;
    end
    for iv=1:nEdge(i)
        if iv>1
            im=i;
        end
        it=iv+cEdge(i);
        itp=it+1;
        if itp>N
            itp=itp-N;
        end
        itm=it-1;
        if itm<1
            itm=itm+N;
        end
        for j=1:nV
            jm=j-1;
            if jm<1
                jm=jm+nV;
            end
            for jv=1:nEdge(j)
                if jv>1
                    jm=j;
                end
                jt=jv+cEdge(j);
                jtp=jt+1;
                if jtp>N
                    jtp=jtp-N;
                end
                jtm=jt-1;
                if jtm<1
                    jtm=jtm+N;
                end
                      Mk11(it,jt)=integral2(@(s,t) PolyInt11(s,t,k,h(i),h(j),u(itp),u(jtp),b(it),b(jt),L,xv(i),yv(i),xv(j),yv(j),CosEdgeAngle(i),SinEdgeAngle(i),CL(i),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(it),b(it),a(jt),b(jt));
                      Mk12(it,jt)=integral2(@(s,t) PolyInt12(s,t,k,h(i),h(j),u(itp),u(jt),b(it),b(jt),L,xv(i),yv(i),xv(j),yv(j),CosEdgeAngle(i),SinEdgeAngle(i),CL(i),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(it),b(it),a(jt),b(jt));
                      Mk21(it,jt)=integral2(@(s,t) PolyInt21(s,t,k,h(i),h(j),u(it),u(jtp),b(it),b(jt),L,xv(i),yv(i),xv(j),yv(j),CosEdgeAngle(i),SinEdgeAngle(i),CL(i),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(it),b(it),a(jt),b(jt));
                      Mk22(it,jt)=integral2(@(s,t) PolyInt22(s,t,k,h(i),h(j),u(it),u(jt),b(it),b(jt),L,xv(i),yv(i),xv(j),yv(j),CosEdgeAngle(i),SinEdgeAngle(i),CL(i),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(it),b(it),a(jt),b(jt));
            
            end
        end
    end
end
  
 Mk=(-1i/4).*(Mk11([end 1:end-1],[end 1:end-1])+Mk12([end 1: end-1], 1:end) + Mk21(1:end,[end 1: end-1])+ Mk22(1:end,1:end));
 M=Mk+(-1i/4)*G0;
 sigma =M\f;            
 
sol=0; %pick from sol=0 to solve at 1 point or sol=1 to solve on many points of a grid

if sol==0
    
 x1 = 0.5; % x cartesian coordinates
 x2 = 0.5;
 ExtInt=zeros(1,N);

for i=1:nV
    im=i-1;
    if im<1
        im=im+nV;
    end
    for iv=1:nEdge(i)
        if iv>1
            im=i;
        end
        it=iv+cEdge(i);
        itm=it-1;
        if itm<1
            itm=itm+N;
        end
        
        ExtInt(it) = (-1i/4)*(integral(@(t) HelmLinSolP(t,u(it),xv(im),yv(im),CL(im),CosEdgeAngle(im),SinEdgeAngle(im),k,h(i),x1,x2,b(itm),L),a(itm),b(itm)) + integral(@(t) HelmLinSolM(t,u(it),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),x1,x2,b(it),L),a(it),b(it)));
     end
 end
 ux = ExtInt*sigma;
 ue = cos(k*x1); % Exact solution
 Eabs = abs(ue-ux); %Absolute Error
 Erel= Eabs/abs(ue) % Relative Error       

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
                     if itp>N
                         itp=itp-N;
                     end
%                      if v1(j)<0 || v2(m)<0
%                        ExtIntP(it) = integral(@(t) HelmLinSolP(t,u(itp),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),v1(j),v2(m),b(it),L),a(it),b(it));
%                        ExtIntM(it) = integral(@(t) HelmLinSolM(t,u(it),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),v1(j),v2(m),b(it),L),a(it),b(it));%  
%                        else
%                          ExtIntP(it)=0;
%                          ExtIntM(it)=0;
%                      end
                       ExtIntP(it) = integral(@(t) HelmLinSolP(t,u(itp),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),v1(j),v2(m),b(it),L),a(it),b(it));
                       ExtIntM(it) = integral(@(t) HelmLinSolM(t,u(it),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),v1(j),v2(m),b(it),L),a(it),b(it));
                       ExtInt= (-1i/4).*(ExtIntP([end 1: end-1]) + ExtIntM(1:end)); %%% Calculate V
                 end
             end
         
             u(j,m)= ExtInt*sigma;
         
          end
      end
   
     mesh(V2,V1,real(u))
     colorbar
     grid off
     title('Numerical solution at k=1 (real part)')
     set(gca,'fontsize',18)
     view(3) 
end

toc










