clear 
tic 
N=16; %no. of elements 
L = 2*pi; %arclength 
s = linspace(0,L,N+1); %element end-points 
a=s(1:N); %integration range
b=s(2:N+1);
v = (cos(a)-cos(b))'; % vector of f(x)
k=1;

for i=1:N
    for j=1:N 
       if i==j 
           A(i,j)=(-1i/4)*(integral2(@(s,si)GkIntMod(s,si,k),a(j),b(j),a(i),@(s)s)+integral2(@(s,si)GkIntMod(s,si,k),a(j),b(j),@(s)s,b(i))+(2i/pi)*(integral(@(si)G0Int(a(j),b(j),si),a(i),b(i)))); %integral on g0int
       else
           A(i,j)=(-1i/4)*integral2(@(s,si)GkInt(s,si,k),a(j),b(j),a(i),b(i)); 
       end
    end
end
sigma =A\v; %finding sigma

sol=0;  % choose 0 for single point, 1 for mesh plot over multiple points

if sol==0

    x1 = 0.5;% x cartesian coordinates
    x2 = 0.5; 


    f = @(s) -(1i/4)*besselh(0,1,k*sqrt((x1-cos(s)).^2 + (x2 - sin(s)).^2)) ; %greens function 
    for i =1:N
        ExtInt(i) = integral(f,a(i),b(i)); %u(x) solution
    end
 
    u = ExtInt*sigma;

    r = sqrt(x1^2 +x2^2); %converting x1, x2 into polar coordinates
    theta = atan2(x2,x1); 

    ue = besselj(1,k*r)*sin(theta)/besselj(1,k); % Exact solution 
    
    Eabs =abs(ue-u); %Absolute Error
    Erel=Eabs/abs(ue); % Relative Error
 toc 
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
        
            f = @(s) -(1i/4)*besselh(0,1,k*sqrt((x1-cos(s)).^2 + (x2 - sin(s)).^2)) ; %greens function 

            for  m=1:N
                 ExtInt(m) = integral(f,a(m),b(m)); %u(x) solution 
            end
            u(i,j) = ExtInt*sigma;
        
            ue(i,j) = besselj(1,k*r(i))*sin(theta(j))/besselj(1,k);
        end
    end
    Eabs =abs(ue-u); %Absolute Error  
    Erel=mean(mean(Eabs))/mean(mean(abs(ue))) % Relative Error

    theta =theta';
 
    V1=cos(theta)*r;
    V2=sin(theta)*r;

    subplot(1,2,1)
    surf(V1',V2',real(u))
    title('Numerical solution (real part)')
    colorbar
    view(2)
    axis equal
    subplot(1,2,2)
    surf(V1',V2',real(ue))
    title('Exact solution (real part)')
    colorbar
    view(2)
    axis equal
end
