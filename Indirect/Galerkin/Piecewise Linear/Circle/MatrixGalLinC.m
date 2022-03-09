    M011=zeros(N,N); M012=zeros(N,N);M021=zeros(N,N); M022=zeros(N,N);
    dx=h;
    for i=1:N
       for j=1:N
            jm=j-1;
            if jm<1
                jm=jm+N;
            end
            im=i-1;
            if im<1
                im=im+N;
            end
            jp=j+1;
            if jp>N
                jp=jp-N;
            end
            ip=i+1;
            if ip>N
                ip=ip-N;
            end
            %%%%%%%%%%% INTEGRAL OF ln(s-t)(s-si+h)/h bj(y) IS SPLIT INTO  ((s-si)ln(s-t)/h + ln(s-t))*bj   
             
            if i~=j
                 M011(i,j)= integral2(@(s,t)G0Num11(t,s,u(ip),b(i),u(jp),b(j),dx,L),a(i),b(i),a(j),b(j));
                 M012(i,j)= integral2(@(s,t)G0Num12(t,s,u(ip),b(i),u(j),b(j),dx,L),a(i),b(i),a(j),b(j)); 
                 M021(i,j)= integral2(@(s,t)G0Num21(t,s,u(i),b(i),u(jp),b(j),dx,L),a(i),b(i),a(j),b(j)); 
                 M022(i,j)= integral2(@(s,t)G0Num22(t,s,u(i),b(i),u(j),b(j),dx,L),a(i),b(i),a(j),b(j)); 
%                  
%                  M011(i,j)= integral2(@(s,t)G0Num11(t,s,u(i),b(im),u(j),b(jm),dx,L),a(im),b(im),a(jm),b(jm));
%                  M012(i,j)= integral2(@(s,t)G0Num12(t,s,u(i),b(im),u(j),b(j),dx,L),a(im),b(im),a(j),b(j)); 
%                  M021(i,j)= integral2(@(s,t)G0Num21(t,s,u(i),b(i),u(j),b(jm),dx,L),a(i),b(i),a(jm),b(jm)); 
%                  M022(i,j)= integral2(@(s,t)G0Num22(t,s,u(i),b(i),u(j),b(j),dx,L),a(i),b(i),a(j),b(j)); 

             elseif i==1 %Calculates the first diagonal entry, the if statement below then                  
                 M011(i,j)=integral2(@(s,t)G0tNum11(s,t,u(j),b(j),dx,L),a(i),b(i),a(j),b(j))+integral(@(t)G0Exact11(a(i),b(i),t,u(j),b(j),dx,L),a(j),b(j));
                 M012(i,j)=integral2(@(s,t)G0tNum12(s,t,u(j),b(j),dx,L),a(i),b(i),a(j),b(j))+integral(@(t)G0Exact12(a(i),b(i),t,u(i),b(i),u(j),b(j),dx,L),a(j),b(j));                 
                 M021(i,j)=integral2(@(s,t)G0tNum21(t,s,u(j),b(j),dx,L),a(i),b(i),a(j),b(j))+integral(@(t)G0Exact21(a(i),b(i),t,u(i),b(i),u(j),b(j),dx,L),a(j),b(j));
                 M022(i,j)=integral2(@(s,t)G0tNum22(t,s,u(j),b(j),dx,L),a(i),b(i),a(j),b(j))+integral(@(t)G0Exact22(a(i),b(i),t,u(j),b(j),dx,L),a(j),b(j));
                 
%                  M011(i,j)=integral2(@(s,t)G0tNum11(s,t,u(j),b(jm),dx,L),a(im),b(im),a(jm),b(jm))+integral(@(t)G0Exact11(a(im),b(im),t,u(j),b(jm),dx,L),a(jm),b(jm));
%                  M012(i,j)=integral2(@(s,t)G0tNum12(s,t,u(j),b(j),dx,L),a(im),b(im),a(j),b(j))+integral(@(t)G0Exact12(a(im),b(im),t,u(i),b(im),u(j),b(j),dx,L),a(j),b(j));                 
%                  M021(i,j)=integral2(@(s,t)G0tNum21(t,s,u(j),b(jm),dx,L),a(i),b(i),a(jm),b(jm))+integral(@(t)G0Exact21(a(i),b(i),t,u(i),b(i),u(j),b(jm),dx,L),a(jm),b(jm));
%                  M022(i,j)=integral2(@(s,t)G0tNum22(t,s,u(j),b(j),dx,L),a(i),b(i),a(j),b(j))+integral(@(t)G0Exact22(a(i),b(i),t,u(j),b(j),dx,L),a(j),b(j));
             end

             if i>1
                 M011(i,i)=M011(1,1);
                 M012(i,i)=M012(1,1);
                 M021(i,i)=M021(1,1);
                 M022(i,i)=M022(1,1);
             end
        end
    end
    M0=(2i/pi).*(M011([end 1:end-1],[end 1:end-1])+M012([end 1: end-1], 1:end) + M021(1:end,[end 1: end-1])+ M022(1:end,1:end));
%     M0=(2i/pi).*(M011+M012+M021+M022);