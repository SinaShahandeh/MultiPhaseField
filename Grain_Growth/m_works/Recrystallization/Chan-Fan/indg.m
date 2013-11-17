% providing symetric or periodic boundary condition
% in any direction if inindex is the grid point index which the field
% value is requested the outindex is the actual point which the field
% should be obtained. gridn is the box size in that direction. For example
% if inindex=gridn+1, the outindex will be first node on the other side of
% the box if we have periodic boundary condition and it will be gridn-1 if
% we want symmetric boundary condition

function [outindex]=indg(inindex,gridn)

%% periodic boundary condition:

outindex=inindex;
if inindex<=0;
    outindex=inindex+gridn;
end
if inindex>gridn
    outindex=inindex-gridn;
end

%% symmetric boundary condition:

% outindex=inindex;
% if inindex<=0;
%     outindex=-inindex+1;
% end
% if inindex>gridn
%     outindex=2*gridn-inindex;
% end
% 


