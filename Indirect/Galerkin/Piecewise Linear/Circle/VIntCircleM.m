function [IntGrand ] = VIntCircleM(t,tj,h)

IntGrand = sin(t).*((tj-t+h)/h); %greens function 

end