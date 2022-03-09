function [Int ] = FCosP(x,si,b,k,h,L)

if si==0 && b==L
    si=L;
end

Int = cos(k.*x).*((x-si+h)/h);
% Int = sin(k.*x*pi).*((x-si+h)/h); %works better

end 