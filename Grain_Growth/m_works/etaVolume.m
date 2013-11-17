% Integrating the eta over the domain to measue each phase amount

function [etaVol]=etaVolume(eta,delx,nboxsize,mboxsize)


 simple summation of points
 for p=1:size(eta,3)
 etaVol(p)=sum(sum(sum(delx^2*eta(:,:,p))));
 end
% more complex spline interpolation and numerical integration

% x=linspace(0,nboxsize*delx,nboxsize);
% y=linspace(0,mboxsize*delx,mboxsize);
% [X,Y]=meshgrid(x,y);
% 
% for p=1:size(eta,3)
%     etaVol(p)=dblquad(@(xi,yi) interp2(X,Y,eta(:,:,p),xi,yi,'cubic'),0,(nboxsize)*delx,0,(mboxsize)*delx,1e-3);
% end


    