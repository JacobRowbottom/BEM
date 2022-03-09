function [IntGrand ] = VIntCircleP(t,tj,b,h,L)
if tj==0 && b==L
    tj=L;
end

IntGrand = sin(t).*((t-tj+h)/h); %greens function 

end