function [IntGrand ] = G0IntPolyNum(t,P,Q)

IntGrand = log(sqrt(t.^2+P.*t+Q)); %greens function 

end