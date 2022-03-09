clear 

SL=zeros(1,4);
nEdge=zeros(1,4);
stepsize=zeros(1,4);
CosEdgeAngle=zeros(1,4);
SinEdgeAngle=zeros(1,4);
CL=zeros(1,4);
cEdge=zeros(1,4);

 %coordinates of unit square
xv = [0 1 1 0]; 
 yv = [0 0 1 1];

%coordinates of L-shape
% xv = [0 1 1 0.5 0.5 0]; 
% yv = [0 0 1 1 0.5 0.5];

nV=length(xv); %no. of vertices of polygons 
NA=8;%number of elements
k=1;

nGaussI=7;
[xg, wg]=grule(nGaussI);

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
Ds = L/NA; %average length of element

for j=1:nV
    nEdge(j) = round(SL(j)/Ds); %Number of elements on each side
    cEdge(j)=sum(nEdge(1:j-1));   % cummaltive number of elements on edge
    stepsize(j)=SL(j)/nEdge(j); %stepsize on each side
    
    for je=1:nEdge(j)
        
        a(je+cEdge(j))=CL(j)+(je-1)*stepsize(j);
        b(je+cEdge(j))=CL(j)+je*stepsize(j);
        
        xa(je+cEdge(j))=xv(j)+stepsize(j)*CosEdgeAngle(j)*(je-1);
        ya(je+cEdge(j))=yv(j)+stepsize(j)*SinEdgeAngle(j)*(je-1); %coordinates 

        xb(je+cEdge(j))=xv(j)+stepsize(j)*CosEdgeAngle(j)*(je);
        yb(je+cEdge(j))=yv(j)+stepsize(j)*SinEdgeAngle(j)*(je);

    end
end

N=sum(nEdge);
for i=1:nV

    for iv=1:nEdge(i)
        it=iv+cEdge(i);
        si=0.5*(b(it)-a(it))*xg'+0.5*(b(it)+a(it));
        wi=0.5*(b(it)-a(it))*wg;
        
%         if i==3                              % vector of L shape
%             v(it,1)=(cos((a(it)-2)*pi)-cos((b(it)-2)*pi))/pi;
%         else 
%             v(it,1)=0;
%         end 
%         
        
    if i==1                                % vector of unit square 
       v(it,1) = (1/k)*(sin(k*b(it))-sin(k*a(it))); 
    elseif i==2
       v(it,1)=cos(k)*(b(it)-a(it)); 
    elseif i==3
        v(it,1) = (1/k)*(sin(k*(3-a(it)))-sin(k*(3-b(it))));
    else
        v(it,1)=(b(it)-a(it)); 
    end
        for j=1:nV
            for jv=1:nEdge(j)
                jt=jv+cEdge(j);
                if it~=jt              
                   M(it,jt)=(-1i/4).*wi*GkIntPolyGal(si,k,xv(i),yv(i),CosEdgeAngle(i),SinEdgeAngle(i),CL(i),stepsize(j),xa(jt),xb(jt),ya(jt),yb(jt),a(jt),b(jt)); % matrix of integrations                
                else                   
                   M(it,jt)= (-1i/4)*(wi*GkIntPolyModD(si,k,xv(i),yv(i),CosEdgeAngle(i),SinEdgeAngle(i),CL(i),stepsize(j),xa(jt),xb(jt),ya(jt),yb(jt),a(jt),b(jt))+(2i/pi)*wi*G0IntPolyD(a(jt),b(jt),si));
                end            
            end 
        end
    end
end

sigma =M\v; %finding sigma

sol=0; %pick from sol=0 to solve at 1 point or sol=1 to solve on many points of a grid

if sol==0
    
x1 = 0.5; % x cartesian coordinates
x2 = 0.5;
ExtInt=zeros(1,N);

 for i =1:nV
    for iv=1:nEdge(i)
            it=iv+cEdge(i);
            ExtInt(it) =integral(@(s)HelmholtzInteriorSol(s,xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,x1,x2),a(it),b(it)); %u(x) solution
    end
end
u = ExtInt*sigma

ue = cos(k*x1); % Exact solution 

Eabs = abs(ue-u); %Absolute Error
Erel= Eabs/abs(ue) % Relative Error

elseif sol==1
    
    nplot=50;
    v1=linspace(-1,1,nplot);
    v2=linspace(-1,1,nplot);
   
   % v1=linspace(0,1,nplot);
   % v2=linspace(0,1,nplot);

    [V1,V2] =meshgrid(v1,v2);
    
    for j=1:nplot;
        j
         for m=1:nplot;
                
             for i =1:nV
%      
                 for iv=1:nEdge(i)
%             
                     it=iv+cEdge(i);
    %                
                     if v1(j)<0 || v2(m)<0
                          ExtInt(it) =integral(@(s)HelmholtzInteriorSol(s,xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,v1(j),v2(m)),a(it),b(it)); %u(x) solution
                     else
                            ExtInt(it)=0;
%                   
                     end
                 end
             end
%         
             u(j,m)= ExtInt*sigma;
%         
          end
      end
%   
     mesh(V2,V1,real(u))
end