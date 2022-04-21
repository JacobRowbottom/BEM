clear 

N=64; %no. of elements 
L = 2*pi; %arclength 
h=L/N; %stepsize 
t= linspace(0,L,N+1); 
a=t(1:N); %element end-points  and collocation points
b=t(2:N+1);
s= linspace(0,L,N+1); 
v = (sin(a))'; % vector of f(x)
k=1;
for i=1:N
    for j=1:N
        jm=j-1;
        if jm<1
           jm=jm+N;
        end

        if i==j
              M(i,j) =(-1i/4)*(integral(@(t)GkLinDp(t,s(i),s(j),k,h),a(jm),b(jm)) +(2i/pi)*(integral(@(t)G0tInt(t,s(i),s(j),h),a(jm),b(jm)) +G0Int(a(jm),b(jm),s(i))) +integral(@(t)GkLinDm(t,s(i),s(j),k,h),a(j),b(j)) +(2i/pi)*(-integral(@(t)G0tInt(t,s(i),s(j),h),a(j),b(j)) +G0Int(a(j),b(j),s(i)))); %mathematical right way 
        else
             M(i,j) = (-1i/4)*(integral(@(t)GkIntLinP(t,s(i),s(j),k,h),a(jm),b(jm))+integral(@(t)GkIntLin(t,s(i),s(j),k,h),a(j),b(j))); 
        end
    end 
end

sigma =M\v; %finding sigma

sol=0;  % choose 0 for single point, 1 for mesh plot over multiple points

if sol==0

    x1 = 0.5;% x cartesian coordinates
    x2 = 0.5; 
  
    for i =1:N
        im=i-1;
        if im<1
           im=im+N;
        end
     g= @(t) besselh(0,1,k*sqrt((x1-cos(t)).^2 + (x2 - sin(t)).^2)).*((s(i)-t+h)/h);
     f= @(t) besselh(0,1,k*sqrt((x1-cos(t)).^2 +(x2 - sin(t)).^2)).*((t -s(i)+h)/h);
        ExtInt(i) = (-1i/4)*(integral(f,a(im),b(im))+ integral(g,a(i),b(i))); 
     
    end
   u = ExtInt*sigma;

    r = sqrt(x1^2 +x2^2); %converting x1, x2 into polar coordinates
    theta = atan2(x2,x1); 

    ue = besselj(1,k*r)*sin(theta)/besselj(1,k); % Exact solution 
    
    Eabs =abs(ue-u); %Absolute Error
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
                f= @(t) besselh(0,1,k*sqrt((x1-cos(t)).^2 +(x2 - sin(t)).^2)).*((t -s(m)+h)/h);
                g= @(t) besselh(0,1,k*sqrt((x1-cos(t)).^2 + (x2 - sin(t)).^2)).*((s(m)-t+h)/h);
                ExtInt(m) = (-1i/4)*(integral(f,a(lm),b(lm))+ integral(g,a(m),b(m))); %u(x) solution 
            end
            u(i,j) =  ExtInt*sigma;
        
            ue(i,j) = besselj(1,k*r(i))*sin(theta(j))/besselj(1,k);
        end
    end
    Eabs =abs(ue-u); %Absolute Error  
    Erel=mean(mean(Eabs))/mean(mean(abs(ue))) % Relative Error

    theta =theta';
 
    V1=cos(theta)*r;
    V2=sin(theta)*r;

end
