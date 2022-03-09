% function [Int ] = FCosMt(x,si,b,k,h,L,CL)
function [Int ] = FCosMt(x,si,b,k,h,L,CL,xv,CosEdgeAngle)

if si==0 && b==L
    si=L;
end
sx=xv+(x-CL)*CosEdgeAngle; 

% Int = cos(k.*(CL-x+1)).*((si-x+h)/h);
Int = cos(k.*sx).*((si-x+h)/h);
end 