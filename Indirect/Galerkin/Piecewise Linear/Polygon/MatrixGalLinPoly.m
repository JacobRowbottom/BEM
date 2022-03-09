
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
                %1 indicates the basis function (t - tj-1)/h 
%2 indicates the basis function (tj+1 - t)/h 
%in Eqn (2.36) in thesis.
%Solved as in Appendix B
                if it==jt                   
                    G011(it,jt)= integral2(@(s,t) G0tIntP(s,t,u(jtp),b(jt),h(i),h(j),L),a(it),b(it),a(jt),b(jt),'RelTol',1e-3,'AbsTol',1e-6)+ integral(@(t)G0ExactPoly11(a(it),b(it),t,u(itp),u(jtp),b(jt),h(i),h(j),L),a(jt),b(jt),'RelTol',1e-3,'AbsTol',1e-6);
                    G012(it,jt)= integral2(@(s,t) G0tIntM(s,t,u(jt),b(jt),h(i),h(j),L),a(it),b(it),a(jt),b(jt),'RelTol',1e-3,'AbsTol',1e-6)+ integral(@(t)G0ExactPoly12(a(it),b(it),t,u(itp),u(jt),b(jt),h(i),h(j),L),a(jt),b(jt),'RelTol',1e-3,'AbsTol',1e-6);
                    G021(it,jt)= -integral2(@(s,t) G0tIntP(s,t,u(jtp),b(jt),h(i),h(j),L),a(it),b(it),a(jt),b(jt),'RelTol',1e-3,'AbsTol',1e-6)+integral(@(t)G0ExactPoly21(a(it),b(it),t,u(it),u(jtp),b(jt),h(i),h(j),L),a(jt),b(jt),'RelTol',1e-3,'AbsTol',1e-6) ; 
                    G022(it,jt)= -integral2(@(s,t) G0tIntM(s,t,u(jt),b(jt),h(i),h(j),L),a(it),b(it),a(jt),b(jt),'RelTol',1e-3,'AbsTol',1e-6)+ integral(@(t)G0ExactPoly22(a(it),b(it),t,u(it),u(jt),b(jt),h(i),h(j),L),a(jt),b(jt),'RelTol',1e-3,'AbsTol',1e-6); 
                else
                    G011(it,jt)=integral2(@(s,t) G0PolyNum11(s,t,h(i),h(j),u(itp),u(jtp),b(it),b(jt),L,xv(i),yv(i),xv(j),yv(j),CosEdgeAngle(i),SinEdgeAngle(i),CL(i),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(it),b(it),a(jt),b(jt),'RelTol',1e-3,'AbsTol',1e-6);
                    G012(it,jt)=integral2(@(s,t) G0PolyNum12(s,t,h(i),h(j),u(itp),u(jt),b(it),b(jt),L,xv(i),yv(i),xv(j),yv(j),CosEdgeAngle(i),SinEdgeAngle(i),CL(i),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(it),b(it),a(jt),b(jt),'RelTol',1e-3,'AbsTol',1e-6);
                    G021(it,jt)=integral2(@(s,t) G0PolyNum21(s,t,h(i),h(j),u(it),u(jtp),b(it),b(jt),L,xv(i),yv(i),xv(j),yv(j),CosEdgeAngle(i),SinEdgeAngle(i),CL(i),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(it),b(it),a(jt),b(jt),'RelTol',1e-3,'AbsTol',1e-6);
                    G022(it,jt)=integral2(@(s,t) G0PolyNum22(s,t,h(i),h(j),u(it),u(jt),b(it),b(jt),L,xv(i),yv(i),xv(j),yv(j),CosEdgeAngle(i),SinEdgeAngle(i),CL(i),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(it),b(it),a(jt),b(jt),'RelTol',1e-3,'AbsTol',1e-6);
                end
            end
        end
    end
    
end

G0=(2i/pi).*(G011([end 1:end-1],[end 1:end-1])+G012([end 1: end-1], 1:end) + G021(1:end,[end 1: end-1])+ G022(1:end,1:end));