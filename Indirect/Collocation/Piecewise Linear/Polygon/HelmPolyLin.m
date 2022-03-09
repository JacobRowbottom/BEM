clear

SL=zeros(1,4);
nEdge=zeros(1,4);
stepsize=zeros(1,4);
CosEdgeAngle=zeros(1,4);
SinEdgeAngle=zeros(1,4);
CL=zeros(1,4);
cEdge=zeros(1,4);


 xv(1)=0; yv(1)=0; %coordinates of vertices 
 xv(2)=1; yv(2)=0;
 xv(3)=1; yv(3)=1;
 xv(4)=0; yv(4)=1;
% 
% xv(1)=0; yv(1)=0; %coordinates of vertices of L shape
% xv(2)=0; yv(2)=1;
% xv(3)=-1; yv(3)=1;
% xv(4)=-1; yv(4)=-1;
% xv(5)=1; yv(5)=-1;
% xv(6)=1; yv(6)=0;

nV=length(xv); %no. of vertices of polygons 

NA=32; %number of elements
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
    stepsize(j)=SL(j)/nEdge(j); %stepsize on each side
    
    for je=1:nEdge(j)
        
        a(je+cEdge(j))=CL(j)+(je-1)*stepsize(j); %element end points 
        b(je+cEdge(j))=CL(j)+je*stepsize(j);
        
        xa(je+cEdge(j))=xv(j)+stepsize(j)*CosEdgeAngle(j)*(je-1);
        ya(je+cEdge(j))=yv(j)+stepsize(j)*SinEdgeAngle(j)*(je-1); %coordinates of elements

        xb(je+cEdge(j))=xv(j)+stepsize(j)*CosEdgeAngle(j)*(je);
        yb(je+cEdge(j))=yv(j)+stepsize(j)*SinEdgeAngle(j)*(je);

    end
end

s=a; xs=xa; ys=ya;

N=sum(nEdge);

% si= 0.5*(a+b); %collocation points 
% xi= 0.5*(xa+xb); %collocation cartesian coordinates
% yi= 0.5*(ya+yb); 


for i=1:nV
     if i==1
       v(1+cEdge(i):cEdge(i)+nEdge(i),1) = cos(k*(xs(1+cEdge(i):cEdge(i)+nEdge(i))))'; % f(x) of square 
    elseif i==2
       v(1+cEdge(i):cEdge(i)+nEdge(i),1)=cos(k)*ones(nEdge(i),1); 
    elseif i==3
        v(1+cEdge(i):cEdge(i)+nEdge(i),1) = cos(k*(xs(1+cEdge(i):cEdge(i)+nEdge(i))))';
    else
        v(1+cEdge(i):cEdge(i)+nEdge(i),1)=ones(nEdge(i),1); 
    end
%     if i==3                             
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
                 dam= sqrt((xs(it)-xa(jm))^2+(ys(it)-ya(jm))^2);
                 dbm= sqrt((xs(it)-xb(jm))^2+(ys(it)-yb(jm))^2);
                 
                 da= sqrt((xs(it)-xa(jt))^2+(ys(it)-ya(jt))^2);
                 db= sqrt((xs(it)-xb(jt))^2+(ys(it)-yb(jt))^2);
                 
                 hm=abs(b(jm)-a(jm));
                 Pm=((dbm^2-dam^2 -hm^2)/hm)-2*a(jm);
                 Qm=dam^2-((dbm^2-dam^2 -hm^2)*a(jm)/hm)+a(jm)^2;
                 
                 P=((db^2-da^2 -stepsize(j)^2)/stepsize(j))-2*a(jt);
                 Q=da^2-((db^2-da^2 -stepsize(j)^2)*a(jt)/stepsize(j))+a(jt)^2;

%                if it~=jt    
                    M(it,jt)=(-1i/4)*(integral(@(t)GkIntPolyP(t,s(jt),b(jm),Pm,Qm,k,hm,L),a(jm),b(jm))+integral(@(t)GkIntPolyM(t,s(jt),P,Q,k,stepsize(j)),a(jt),b(jt)));
%                    %M(it,jt)=(-1i/4)*integral(@(t)GkIntPolyP(t,s(jt),b(jm),Pm,Qm,k,hm,L),a(jm),b(jm));
%                else
%                   M(it,jt) = (-1i/4)*(integral(@(t)GkPolyDP(t,s(jt),b(jm),Pm,Qm,k,hm,L),a(jm),b(jm)) +(2i/pi)*(integral(@(t)G0tPolyInt(t,s(jt),a(jm),b(jm),Pm,Qm,hm,L),a(jm),b(jm))+ G0IntPolyD(a(jm),b(jm),s(it),L))+ integral(@(t)GkPolyDM(t,s(jt),P,Q,k,stepsize(j)),a(jt),b(jt)) + (2i/pi)*(-integral(@(t)G0tPolyInt(t,s(jt),a(jt),b(jt),P,Q,stepsize(j),L),a(jt),b(jt))+ G0IntPolyD(a(jt),b(jt),s(it),L)));
%                end
%             
            end 
        end
    end
end

sigma =M\v; %finding sigma
    
x1 = 0.5; % x cartesian coordinates
x2 = 0.5;
ExtInt=zeros(1,N);

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
                     
%         ExtInt(i) = (-1i/4)*(integral(f,a(im),b(im))+ integral(g,a(i),b(i))); % right way 
          if iv>1
              
              ExtInt(it) = (-1i/4)*(integral(@(t)HelmLinSolP(t,s(it),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,stepsize(i),x1,x2,b(im),L),a(im),b(im)) + integral(@(t)HelmLinSolM(t,s(it),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,stepsize(i),x1,x2),a(it),b(it))); %u(x) solution
          else
              
              ExtInt(it) = (-1i/4)*(integral(@(t)HelmLinSolP(t,s(it),xv(ii),yv(ii),CL(ii),CosEdgeAngle(ii),SinEdgeAngle(ii),k,hm,x1,x2,b(im),L),a(im),b(im)) + integral(@(t)HelmLinSolM(t,s(it),xv(i),yv(i),CL(i),CosEdgeAngle(i),SinEdgeAngle(i),k,stepsize(i),x1,x2),a(it),b(it))); %u(x) solution
             
          end
              
              
          end
   end
   u = ExtInt*sigma;
   ue = cos(k*x1); % Exact solution 

Eabs = abs(ue-u); %Absolute Error
Erel= Eabs/abs(ue) % Relative Error

