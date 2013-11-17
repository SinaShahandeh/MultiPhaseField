function [outindex]=indg(inindex,gridn)

%% Periodic boundary condition:
 outindex=inindex;
 if inindex<=0;
     outindex=inindex+gridn;
 end
 if inindex>gridn
     outindex=inindex-gridn;
 end

%% Symmetric boundary condition:
%outindex=inindex;
%if inindex<=0;
%    outindex=-inindex+1;
%end
%if inindex>gridn
%    outindex=2*gridn-inindex;
%end

