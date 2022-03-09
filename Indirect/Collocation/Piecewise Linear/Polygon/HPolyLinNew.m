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

NA=32; %number of elements
k=10;   %wavenumber  

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

N=sum(nEdge);

%%% Asin runs script to solve M2 which solves the martix entries for the
%%% ingetration of G0 in the singularity subtraction procedure
%   Asin;
%  Bsin;
for i=1:nV
     if i==1
       v(1+cEdge(i):cEdge(i)+nEdge(i),1) = cos(k*(xs(1+cEdge(i):cEdge(i)+nEdge(i))))'; % boundary conditions for square 
    elseif i==2
       v(1+cEdge(i):cEdge(i)+nEdge(i),1)=cos(k)*ones(nEdge(i),1); 
    elseif i==3
        v(1+cEdge(i):cEdge(i)+nEdge(i),1) = cos(k*(xs(1+cEdge(i):cEdge(i)+nEdge(i))))';
    else
        v(1+cEdge(i):cEdge(i)+nEdge(i),1)=ones(nEdge(i),1); 
     end
%         if i==6
%             v(1+cEdge(i):cEdge(i)+nEdge(i),1)=cos(k*(xs(1+cEdge(i):cEdge(i)+nEdge(i))))'; % vector of L shape
%         else 
%             v(1+cEdge(i):cEdge(i)+nEdge(i),1)=0;
%         end 
     
      for iv=1:nEdge(i)
         it=iv+cEdge(i);
         for j=1:nV
            for jv=1:nEdge(j)
                jt=jv+cEdge(j);
                   jm=jt-1;
                   if jm<1
                      jm=jm+N;
                   end
                   
                   dam= sqrt((xs(it)-xa(jm))^2+(ys(it)-ya(jm))^2); %distance from collaction point and element end points, for previous element
                   dbm= sqrt((xs(it)-xb(jm))^2+(ys(it)-yb(jm))^2);
                   
                   da= sqrt((xs(it)-xa(jt))^2+(ys(it)-ya(jt))^2);%distance from collaction point and element end points, for element
                   db= sqrt((xs(it)-xb(jt))^2+(ys(it)-yb(jt))^2);
                   
                   hm=abs(b(jm)-a(jm)); %stepsize of previous elements 
                   Pm=((dbm^2-dam^2 -hm^2)/hm)-2*a(jm); %Pm,Qm are defined for the distance of collocation point to point t, previous element
                   Qm=dam^2-((dbm^2-dam^2 -hm^2)*a(jm)/hm)+a(jm)^2;
                 
                   P=((db^2-da^2 -h(j)^2)/h(j))-2*a(jt);%P,Q are defined for the distance of collocation point to point t in the element
                   Q=da^2-((db^2-da^2 -h(j)^2)*a(jt)/h(j))+a(jt)^2;
                   %%%%% Mk calculate the integral of (Gk -G0) times basis functions as part of singularity subtraction                  
                   Mk(it,jt)=(-1i/4)*(integral(@(t)GkPolyDP(t,s(jt),b(jm),Pm,Qm,k,hm,L),a(jm),b(jm))+integral(@(t)GkPolyDM(t,s(jt),P,Q,k,h(j)),a(jt),b(jt)));
                   
%                    Mk_P(it,jt)=(-1i/4)*(integral(@(t)GkPolyDP(t,s(jt),b(j),P,Q,k,stepsize(j),L),a(jt),b(jt)));
%                    Mk_M(it,jt)=(-1i/4)*(integral(@(t)GkPolyDM(t,s(jt),P,Q,k,stepsize(j)),a(jt),b(jt)));
            end
        end
    end
end
% MkTest = Mk_P(1:end,[end 1:end-1]) + Mk_M;
M=Mk+(-1i/4)*M2; %%%adds the matrix entries for Mk and M2 

sigma =M\v; %calculates sigma
sol=1; %pick from sol=0 to solve at 1 point or sol=1 to solve on many points of a grid

if sol==0
    
x1 = 0.75; % x cartesian coordinates
x2 = 0.25;
uInter=zeros(1,N);

   for i =1:nV
       ii=i-1;
       if ii<1
          ii=ii+nV;
       end
            
    for iv=1:nEdge(i)
            it=iv+cEdge(i);
            
            im=it-1;
            if im<1
                im=im+N;
            end
            
%             if iv>1
%                 uInter(it) = (-1i/4)*(integral(@(t)HelmLinSolP(t,s(it),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),x1,x2,b(im),L),a(im),b(im)) + integral(@(t)HelmLinSolM(t,s(it),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),x1,x2),a(it),b(it)));
%             else
                uInter(it) = (-1i/4)*(integral(@(t)HelmLinSolP(t,s(it),xv(ii),yv(ii),CL(ii),CosEdgeAngle(ii),SinEdgeAngle(ii),k,hm,x1,x2,b(im),L),a(im),b(im)) + integral(@(t)HelmLinSolM(t,s(it),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),x1,x2),a(it),b(it)));
%             end
    end
   end
   u = uInter*sigma

   ue = cos(k*x1); % Exact solution 

Eabs = abs(ue-u); %Absolute Error
% Erel= Eabs/abs(ue) % Relative Error

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
%                  im = i-1;
%                  if im<1
%                     im=im+nV;
%                  end    
                 for iv=1:nEdge(i)
                     it=iv+cEdge(i);
%                      itm=it-1;
%                      if itm<1
%                          itm=itm+N;
%                      end
%                      uInter(it) = (-1i/4)*(integral(@(t) HelmLinSolP(t,s(it),xv(im),yv(im),CL(im),CosEdgeAngle(im),SinEdgeAngle(im),k,h(im),v1(j),v2(m),b(itm),L),a(itm),b(itm)) + integral(@(t) HelmLinSolM(t,s(it),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),v1(j),v2(m)),a(it),b(it)));
                     uInterP(it) = (-1i/4)*integral(@(t) HelmLinSolP(t,s(it),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),v1(j),v2(m),b(it),L),a(it),b(it));
                     uInterM(it) = (-1i/4)*integral(@(t) HelmLinSolM(t,s(it),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,h(i),v1(j),v2(m)),a(it),b(it));

                 end
             end
             uInterPF = [uInterP(end),uInterP(1:end-1)];
             uInter = uInterP +uInterM;
             u(j,m)= (0.5)*uInter*sigma;
             ue(j,m) = cos(k*v1(j));
          end
      end
   
     mesh(V2,V1,real(u))
%      colorbar
%      grid off
%      title('Numerical solution at k=5 (real part)')
%      set(gca,'fontsize',18)
%      view(3)

     figure
     mesh(V2,V1,real(ue))
end
% save('u512','u');              
                 
                 
                 
                 
                 
                 
                 
                 