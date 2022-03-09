function [Int ] = FCosM(x,si,b,k,h,L)

if si==0 && b==L
    si=L;
end

Int = cos(k.*x).*((si-x+h)/h);
% Int = sin(k.*pi*x).*((si-x+h)/h); %Works better with this one

end 