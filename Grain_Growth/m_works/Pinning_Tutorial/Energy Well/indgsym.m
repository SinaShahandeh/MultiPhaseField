function [outindex]=indgsym(inindex,gridn)

%% Symmetric boundary condition:
outindex=inindex;
if inindex<=0;
    outindex=-inindex+1;
end
if inindex>gridn
    outindex=2*gridn-inindex;
end
