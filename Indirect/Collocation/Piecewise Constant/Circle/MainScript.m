clear 

N=16; %no. of elements 
L = 2*pi; %Perimeter 
s = linspace(0,L,N+1); %arclength
a=s(1:N); %Element end points
b=s(2:N+1);
si= 0.5*(a+b); %collocation points 
f = sin(si)'; % Boundary condition f(x)
k=1;

%Matrix calculations of Eqn (2.30)
for i=1:N
    for j=1:N 
       if i==j %Perform Singularity Subtraction when i=j
           Gk(i,j)=(-1i/4)*(integral(@(s)GkIntMod(s,si(i),k),a(j),b(j))+(2i/pi)*G0Int(a(j),b(j),si(i))); %Singularity Subtraction
       else
           Gk(i,j)=(-1i/4)*integral(@(s)GkInt(s,si(i),k),a(j),b(j));
          
       end
    end
end
sigma =Gk\f; %Calculate Sigma

sol=1;  % choose 0 for single point, 1 for mesh plot over multiple points
if sol==0
   % Cartesian coordinates of interior point
    x1 = 0.5;
    x2 = 0.5; 
    f = @(s) -(1i/4)*besselh(0,1,k*sqrt((x1-cos(s)).^2 + (x2 - sin(s)).^2)) ; %greens function 
    for i =1:N
        Int(i) = integral(f,a(i),b(i)); 
    end
 
    u = Int*sigma; %u(x) solution

    r = sqrt(x1^2 +x2^2); %converting x1, x2 into polar coordinates
    theta = atan2(x2,x1); 

    ue = besselj(1,k*r)*sin(theta)/besselj(1,k); % Exact solution 
       
elseif sol==1
   
    nr=20;
    ntheta=20;
    r=linspace(0,1,nr);
    theta=linspace(-pi,pi,ntheta);

    for i=1:nr;
        i
        for j=1:ntheta;
            x1 = r(i)*cos(theta(j));
            x2 = r(i)*sin(theta(j));
        
            f = @(s) -(1i/4)*besselh(0,1,k*sqrt((x1-cos(s)).^2 + (x2 - sin(s)).^2)) ; %greens function 

            for  m=1:N
                 Int(m) = integral(f,a(m),b(m)); 
            end
            u(i,j) = Int*sigma;%u(x) solution 
        
            ue(i,j) = besselj(1,k*r(i))*sin(theta(j))/besselj(1,k);%Exact solution
        end
    end
    theta =theta';
 
    V1=cos(theta)*r;
    V2=sin(theta)*r;

    figure
%     subplot(1,2,1)
    surf(V2,V1,real(u).')
    title('Numerical solution (real part)')
    colorbar
    view(2)
    axis equal
    figure
%     subplot(1,2,2)
    surf(V2,V1,real(ue).')
    title('Exact solution (real part)')
    colorbar
    view(2)
    axis equal
end

