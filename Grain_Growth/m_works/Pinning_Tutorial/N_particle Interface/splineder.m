% calculates drivative of the spline sp and gives spd another spline
function [spd]=splineder(sp)
% breaking the spline:
[breaks,coefs,l,k,d]=unmkpp(sp);
coefs_d=zeros(size(coefs,1),size(coefs,2)-1);
% because stupid matlab may return 2 degree sometimes for derivative of
% coefs:
coefs=coefs+eps/1e50;
% calculation of derivatives
for ni=1:l
    coefs_d(ni,:)=polyder(coefs(ni,:));
end

% again making spline from derivatives:

spd=mkpp(breaks,coefs_d);
