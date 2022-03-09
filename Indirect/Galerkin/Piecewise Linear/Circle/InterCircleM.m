function [IntGrand ] = InterCircleM(k,x1,x2,t,tj,b,h,L)
if tj==0 && b==L
    tj=L;
end

IntGrand = besselh(0,1,k*sqrt((x1-cos(t)).^2 + (x2 - sin(t)).^2)).*((tj-t+h)/h); %greens function 

end