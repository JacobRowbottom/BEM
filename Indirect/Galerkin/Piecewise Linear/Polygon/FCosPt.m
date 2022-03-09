% function [Int ] = FCosPt(x,si,b,k,h,L,CL)
function [Int ] = FCosPt(x,si,b,k,h,L,CL,xv,CosEdgeAngle)

if si==0 && b==L
    si=L;
end
sx=xv+(x-CL)*CosEdgeAngle; 

% Int = cos(k.*(CL-x+1)).*((x-si+h)/h);
Int = cos(k.*sx).*((x-si+h)/h);
end 