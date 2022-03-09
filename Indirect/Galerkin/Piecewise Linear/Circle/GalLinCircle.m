clear

tic
N=32; %no. of elements 
L = 2*pi; %arclength 
h=L/N; %stepsize 
%s= linspace(0,L,N+1); 
%t= linspace(0,L,N+1); %element end-points 
u= linspace(0,L,N+1); 
a=u(1:N); %integration range
b=u(2:N+1);

k=1;

%%%%%% Calculate the vector f(x) of boundary condition sin(x)
for i =1:N
        im=i-1;
        if im<1
           im=im+N;
        end
        %Integral of Boundary condition f(x)=sin(x) for each basis function 
        f(i) = integral(@(t)VIntCircleP(t,u(i),b(im),h,L),a(im),b(im))+ integral(@(t)VIntCircleM(t,u(i),h),a(i),b(i));
                       
end

%%%%%%%%% INTEGRATE ALL ELEMENTS USING SINGULARITY SUBTRACTION : (Gk -G0)bibj dxdy + iint(G0 bibj) dxdy
%Calculates entries via Appendix B
MatrixGalLinC; %Calculates M0=iint(G0 bibj) dxdy 


%1 indicates the basis function (t - tj-1)/h 
%2 indicates the basis function (tj+1 - t)/h 
%in Eqn (2.36) in thesis.
Mk11=zeros(N,N); %Matrix with basis functions 
Mk12=zeros(N,N);
Mk21=zeros(N,N); 
Mk22=zeros(N,N);
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
        %%%%% Calculates iint (Gk -G0)bibj dxdy
        Mk11(i,j) = integral2(@(s,t) Int11(s,t,u(i),u(j),b(i),b(j),k,h,L),a(i),b(i),a(j),b(j));
        Mk12(i,j) = integral2(@(s,t) Int12(s,t,u(i),u(j),b(i),b(j),k,h,L),a(i),b(i),a(j),b(j));
        Mk21(i,j) = integral2(@(s,t) Int21(s,t,u(i),u(j),b(i),b(j),k,h,L),a(i),b(i),a(j),b(j));
        Mk22(i,j) = integral2(@(s,t) Int22(s,t,u(i),u(j),b(i),b(j),k,h,L),a(i),b(i),a(j),b(j));            
    end 
end
%Add the matrices entried to the correct element entry
    Mk=(-1i/4).*(Mk11([end 1:end-1],[end 1:end-1])+Mk12([end 1: end-1], 1:end) + Mk21(1:end,[end 1: end-1])+ Mk22(1:end,1:end));
    M=Mk+(-1i/4)*M0;
sigma =M\f'; %finding sigma

sol=0;  % choose 0 for single point, 1 for mesh plot over multiple points

if sol==0
    x1 = 0.25;% x cartesian coordinates
    x2 = 0.5; 
    for i =1:N
        im=i-1;
        if im<1
           im=im+N;
        end

        ExtInt(i) = (-1i/4)*(integral(@(t)InterCircleP(k,x1,x2,t,u(i),b(im),h,L),a(im),b(im)) + integral(@(t)InterCircleM(k,x1,x2,t,u(i),b(i),h,L),a(i),b(i))); % right way 
    end
   ua = ExtInt*sigma % approx solution 

    r = sqrt(x1^2 +x2^2); %converting x1, x2 into polar coordinates
    theta = atan2(x2,x1); 

    ue = besselj(1,k*r)*sin(theta)/besselj(1,k); % Exact solution 
%     ue = besselj(3,k*r)*cos(3*theta)/besselj(3,k); 
    
    Eabs =abs(ue-ua) %Absolute Error
    Erel=Eabs/abs(ue) % Relative Error
 
elseif sol==1
   
    nr=20;
    ntheta=40;
    r=linspace(0,1,nr);
    theta=linspace(-pi,pi,ntheta);

    for i=1:nr;
        i
        for j=1:ntheta;
            x1 = r(i)*cos(theta(j));
            x2 = r(i)*sin(theta(j));

            for  m=1:N
           lm=m-1;
        if lm<1
           lm=lm+N;
        end
                f= @(t) besselh(0,1,k*sqrt((x1-cos(t)).^2 +(x2 - sin(t)).^2)).*((t -u(m)+h)/h);
                g= @(t) besselh(0,1,k*sqrt((x1-cos(t)).^2 + (x2 - sin(t)).^2)).*((u(m)-t+h)/h);
                ExtInt(m) = (-1i/4)*(integral(f,a(lm),b(lm))+ integral(g,a(m),b(m))); %u(x) solution 
            end
            uNum(i,j) =  ExtInt*sigma;
        
            ue(i,j) = besselj(1,k*r(i))*sin(theta(j))/besselj(1,k);
        end
    end
%     Eabs =abs(ue-u); %Absolute Error  
%     Erel=mean(mean(Eabs))/mean(mean(abs(ue))) % Relative Error

    theta =theta';
 
    V1=cos(theta)*r;
    V2=sin(theta)*r;

    subplot(1,2,1)
    surf(V2,V1,real(uNum).')
    title('Numerical solution (real part)')
    colorbar
    view(2)
    axis equal
    subplot(1,2,2)
    surf(V2,V1,real(ue).')
    title('Exact solution (real part)')
    colorbar
    view(2)
    axis equal
end
toc






