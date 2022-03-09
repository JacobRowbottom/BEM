function [Int ] = G0PolyNum11(s,t,hi,hj,si,tj,bi,bj,L,sxv,syv,txv,tyv,CosEdgeAngleS,SinEdgeAngleS,CLS,CosEdgeAngleT,SinEdgeAngleT,CLT)
if si==0 && bi==L
    si=L;
end
if tj==0 && bj==L
    tj=L;
end

sx=sxv+(s-CLS)*CosEdgeAngleS; %cartesian coordinates in terms of arclength
sy=syv+(s-CLS)*SinEdgeAngleS;

tx=txv+(t-CLT)*CosEdgeAngleT; %cartesian coordinates in terms of arclength
ty=tyv+(t-CLT)*SinEdgeAngleT;

D=sqrt((sx-tx).^2+(sy-ty).^2);
Ker =log(D); 

Int=Ker.*((s-si+hi)/hi).*((t-tj+hj)/hj);

% s
% sx
% sy
% sxv
% syv
% 
% t
% tx
% ty
% 
% si
% tj
% pause
% 


end