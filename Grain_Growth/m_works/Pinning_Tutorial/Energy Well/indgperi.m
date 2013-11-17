function [outindex]=indgperi(inindex,gridn)

%% Periodic boundary condition:
outindex=inindex;
if inindex<=0;
    outindex=inindex+gridn;
end
if inindex>gridn
    outindex=inindex-gridn;
end


