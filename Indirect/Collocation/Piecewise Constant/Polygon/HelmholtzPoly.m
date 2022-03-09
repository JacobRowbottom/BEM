
clear 

SL=zeros(1,4);
nEdge=zeros(1,4);
h=zeros(1,4);
CosEdgeAngle=zeros(1,4);
SinEdgeAngle=zeros(1,4);
CL=zeros(1,4);
cEdge=zeros(1,4);
%coordinates of vertices 
%Unit Square
xv=[0 1 1 0];
yv=[0 0 1 1];
%L-shape
%xv = [0 1 1 0.5 0.5 0];
% yv = [0 0 1 1 0.5 0.5];
NVert=length(xv); %no. of vertices of polygons 

N=512; %number of elements
k=1; %Wavenumber

for j=1:NVert
    jp = j+1;
    if (jp > NVert)
        jp=jp-NVert;
    end
    SL(j)=sqrt((xv(j)-xv(jp))^2+(yv(j)-yv(jp))^2); %side lengths 
    CL(j)=sum(SL(1:j-1)); %cummalative length 
    CosEdgeAngle(j)=(xv(jp)-xv(j))/SL(j);
    SinEdgeAngle(j)=(yv(jp)-yv(j))/SL(j);
    
end

L=sum(SL); %perimeter of polygon
Ds = L/N; %average length of element

for j=1:NVert
    nEdge(j) = round(SL(j)/Ds); %Number of elements on each side
    cEdge(j)=sum(nEdge(1:j-1));   % cummaltive number of elements on edge
    h(j)=SL(j)/nEdge(j); %stepsize on each side
    
    for je=1:nEdge(j)
        
        a(je+cEdge(j))=CL(j)+(je-1)*h(j);
        b(je+cEdge(j))=CL(j)+je*h(j);
        
        ax(je+cEdge(j))=xv(j)+h(j)*CosEdgeAngle(j)*(je-1);
        ay(je+cEdge(j))=yv(j)+h(j)*SinEdgeAngle(j)*(je-1); %coordinates 

        bx(je+cEdge(j))=xv(j)+h(j)*CosEdgeAngle(j)*(je);
        by(je+cEdge(j))=yv(j)+h(j)*SinEdgeAngle(j)*(je);

    end
end

N=sum(nEdge);

CPi= 0.5*(a+b); %collocation points 
xi= 0.5*(ax+bx); %collocation cartesian coordinates
yi= 0.5*(ay+by); 

%Edit if domain is L-shape. 
for i=1:NVert
    for ie=1:nEdge(i)
       it = ie +cEdge(i);
     %%% f(x) Dirichlet BC
       if i==1
           f(it,1) = cos(k*(xi(it))); %  square 
       elseif i==2
           f(it,1)=cos(k); 
       elseif i==3
           f(it,1) = cos(k*(xi(it)));
       else
           f(it,1)=1; 
       end
     % f(x) of L shape
%        if i==6                             
%            f(it,1)=cos(k*(xi(it)));
%        else 
%            f(it,1)=0;
%        end 

       r = sqrt((bx(it)-ax(it))^2 + (by(it)-ay(it))^2); 
       txq(it) = (bx(it) - ax(it))/r;
       tyq(it) = (by(it) - ay(it))/r;  
    end
end
nxq = -tyq; %% cartesian coordinates of unit normal vector of q. 
nyq = txq;

VSLP = zeros(N,N);  %%%% Matrix calculating for the Dirichlet problem
for i=1:NVert
    for iv=1:nEdge(i)
        it=iv+cEdge(i);
        for j=1:NVert
            for jv=1:nEdge(j)
                jt=jv+cEdge(j);
                if it==jt    
                    VSLP(it,jt)=(-1i/4)*(integral(@(q)VSLPIntPolyDiagonal(q,k,xi(it),yi(it),xv(j),yv(j),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(jt),b(jt)) + (2i/pi)*G0IntPolyD(a(jt),b(jt),CPi(it)));
                else
                    VSLP(it,jt)=(-1i/4)*integral(@(q)VSLPIntPoly(q,k,xi(it),yi(it),xv(j),yv(j),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(jt),b(jt));  
                end
            end
        end
    end
end

sigma =VSLP\f; % sigma for Dirichlet 
sol=0; %pick from sol=0 to solve at 1 point or sol=1 to solve on many points of a grid

if sol==0
    
x1 = 0.5; % x cartesian coordinates
x2 = 0.5;
uInter=zeros(1,N);

 for i =1:NVert
     for iv=1:nEdge(i)
            it=iv+cEdge(i);
            uInter(it) =integral(@(s)HelmholtzInteriorSol(s,xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,x1,x2),a(it),b(it)); %u(x) solution
    end
end
u = uInter*sigma;
uBound =VSLP*sigma;

ue = cos(k*x1); % Exact solution dirichlet
% uBoundE = (1/k)*(cot(k)*cos(k.*xi) + sin(k.*xi)).';% Exact solution Neumann

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
             for i =1:NVert
                 for iv=1:nEdge(i)
                     it=iv+cEdge(i);
                     %Solution for L-shape
%                      if v1(j)<0.5 && v2(m)>0.5
%                          uInter(it)=0;                    
%                      else
%                          uInter(it) =integral(@(s)HelmholtzInteriorSol(s,xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,v1(j),v2(m)),a(it),b(it)); %u(x) solution
%                      end
                     
                     uInter(it) =integral(@(s)HelmholtzInteriorSol(s,xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,v1(j),v2(m)),a(it),b(it)); %u(x) solution
                   
                 end
             end
%              ue(j,m) = cos(k*v1(j));
             u(j,m)= uInter*sigma;
         
          end
    end
end