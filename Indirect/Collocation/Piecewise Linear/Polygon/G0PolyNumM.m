function [Int ] = G0PolyNumM(t,h,tj,bj,L,txv,tyv,CosEdgeAngleT,SinEdgeAngleT,CLT)

if tj==0 && bj==L
    tj=L;
end

tx=txv+(t-CLT)*CosEdgeAngleT; %cartesian coordinates in terms of arclength
ty=tyv+(t-CLT)*SinEdgeAngleT;

D=sqrt((sx-tx).^2+(sy-ty).^2);
Ker =log(D); 

Int=Ker.*((t-tj+h)/h);

end