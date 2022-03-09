%Calculating the Matrix M2 
 for i=1:nV
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
                 
                 P=((db^2-da^2 -h(j)^2)/h(j))-2*a(jt);
                 Q=da^2-((db^2-da^2 -h(j)^2)*a(jt)/h(j))+a(jt)^2;
               
                 
                %%%% We split the integral of G0(t-tj +h)/h into
                %%%% integral of G0(t-tj) [Solved Numerically] +G0 [solved
                %%%% primarily analytically]
                   M0t(it,jt) = (integral(@(t)G0tPolyInt(t,s(jt),a(jm),b(jm),Pm,Qm,hm,L),a(jm),b(jm))-integral(@(t)G0tPolyInt(t,s(jt),a(jt),b(jt),P,Q,h(j),L),a(jt),b(jt))); 
                 
                if it==jt
                    M0(it,jt) = G0IntPolyD(a(jm),b(jm),s(it),L)+G0IntPolyD(a(jt),b(jt),s(it),L);
                elseif s(it)==b(jt)
                    M0(it,jt) = integral(@(t)G0IntPolyNum(t,Pm,Qm),a(jm),b(jm))+ G0IntPolyD(a(jt),b(jt),s(it),L);
                elseif s(it)==0 && b(jt)==L
                    M0(it,jt) = integral(@(t)G0IntPolyNum(t,Pm,Qm),a(jm),b(jm))+ G0IntPolyD(a(jt),b(jt),s(it),L);
                elseif s(it)==a(jm) 
                    M0(it,jt) = G0IntPolyD(a(jm),b(jm),s(it),L)+integral(@(t)G0IntPolyNum(t,P,Q),a(jt),b(jt));
                else 
                    M0(it,jt) =  integral(@(t)G0IntPolyNum(t,Pm,Qm),a(jm),b(jm))+integral(@(t)G0IntPolyNum(t,P,Q),a(jt),b(jt));
                end
            
            end
        end
    end
end
M2=(2i/pi)*(M0+M0t);