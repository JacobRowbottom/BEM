clear
tic

%% Prelimnaries and set up 
SL=zeros(1,4);
nEdge=zeros(1,4);
h=zeros(1,4);
CosEdgeAngle=zeros(1,4);
SinEdgeAngle=zeros(1,4);
CL=zeros(1,4);
cEdge=zeros(1,4);

 % Unit square
 xv=[0 1 1 0];
 yv=[0 0 1 1];
% L- shape
% xv = [0 1 1 0.5 0.5 0]; 
% yv = [0 0 1 1 0.5 0.5];

M=256; %number of elements
k=1;   %wavenumber  
PSx = [0.5]; 
PSy = [0.5]; 

NVert=length(xv); %no. of vertices of polygons 
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
dx = L/M; %average length of element/ 

F_BC = zeros(M,1);
txq = zeros(M,1);
tyq = zeros(M,1);
for j=1:NVert
    nEdge(j) = round(SL(j)/dx); %Number of elements on each side
    cEdge(j)=sum(nEdge(1:j-1));   % cummaltive number of elements on edge
    h(j)=SL(j)/nEdge(j); %stepsize on each side
    
   for je=1:nEdge(j)
        a(je+cEdge(j))=CL(j)+(je-1)*h(j); %element end points 
        b(je+cEdge(j))=CL(j)+je*h(j);
       
        ax(je+cEdge(j))=xv(j)+h(j)*CosEdgeAngle(j)*(je-1);
        ay(je+cEdge(j))=yv(j)+h(j)*SinEdgeAngle(j)*(je-1); %coordinates of elements

        bx(je+cEdge(j))=xv(j)+h(j)*CosEdgeAngle(j)*(je);
        by(je+cEdge(j))=yv(j)+h(j)*SinEdgeAngle(j)*(je);
   end
   
   for je=1:nEdge(j)
       jt = je +cEdge(j);     
       if j==NVert %% f(x)= du/dn, 
           F_BC(jt,1) = 1; 
       else 
           F_BC(jt,1) = 0; 
       end 
       r = sqrt((bx(jt)-ax(jt))^2 + (by(jt)-ay(jt))^2); 
       txq(jt) = (bx(jt) - ax(jt))/r;
       tyq(jt) = (by(jt) - ay(jt))/r;  
   end
end
nx = -tyq; %% cartesian coordinates of unit normal vector of q. 
ny = txq;

CPi= 0.5*(a+b); %collocation points 
xi= 0.5*(ax+bx); %collocation cartesian coordinates
yi= 0.5*(ay+by); 

%% Boundary Solution Calculation
 VSLP = zeros(M,M);  
 KDLP = zeros(M,M);

     for i=1:NVert
         for iv=1:nEdge(i)
             it=iv+cEdge(i);
             for j=1:NVert
                 for jv=1:nEdge(j)
                     jt=jv+cEdge(j);
                        if it==jt    
                            VSLP(it,jt)=(-1i/4)*(integral(@(q)VSLPIntPolyDiagonal(q,k,xi(it),yi(it),xv(j),yv(j),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(jt),b(jt)) + (2i/pi)*G0IntPolyD(a(jt),b(jt),CPi(it)));
%                           KDLP(it,jt) = 0;
                        else
                            VSLP(it,jt)=(-1i/4)*integral(@(q)VSLPIntPoly(q,k,xi(it),yi(it),xv(j),yv(j),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(jt),b(jt));  
                            KDLP(it,jt)=(-1i/4)*integral(@(q)KIntPoly(q,k,xi(it),yi(it),xv(j),yv(j),CosEdgeAngle(j),SinEdgeAngle(j),CL(j),nx(jt),ny(jt)),a(jt),b(jt));      
                        end
                 end
             end
         end
     end
     KLHS = ((0.5)*eye(M) + KDLP); %LHS of Helm eqs 
     VRHS = VSLP*F_BC; %RHS of Helm eqs
     uBound = KLHS\VRHS; %Boundary solution

     %Exact Boundary solution
% uExactBound = (1/k)*(cot(k)*cos(k.*xi) + sin(k.*xi)).';

%% INTERIOR SOLUTION CALCULATION 
stage=2
InteriorNum; 
%% POINT SOURCE CALCULATION 
stage =3
PS_BC = zeros(M,1);
for i=1:NVert
    for iv=1:nEdge(i)
        it=iv+cEdge(i);
        PS_BC(it) = (1i/4)*dGkdn(k,xi(it),yi(it),nx(it),ny(it),PSx,PSy);
    end
end

VLHS= VSLP*PS_BC;
PSBound = KLHS\VLHS;

%%%Exact Sol of point source
PointSourceNum;
PointSourceExact; 


toc
