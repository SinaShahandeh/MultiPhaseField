% Integrating the eta over the domain to measue each phase amount

function [etaVol]=etaVolume(etai,delx,nboxsize,mboxsize,accuracy)

% simple summation of points
% for p=1:size(eta,3)
% etaVol(p)=sum(sum(sum(delx^2*eta(:,:,p))));
% end
% more complex spline interpolation and numerical integration

% for 2D domain
if size(etai,2)>1
    if strcmp(accuracy,'high')

        x=linspace(0,nboxsize*delx,nboxsize);
        y=linspace(0,mboxsize*delx,mboxsize);
        [X,Y]=meshgrid(x,y);
        etaVol=dblquad(@(xi,yi) interp2(X,Y,etai,xi,yi,'cubic'),0,(nboxsize)*delx,0,(mboxsize)*delx,1e-6);

    else
        etaVol=sum(sum(etai))*delx*delx;
    end
end

% for 1D domain
if size(etai,2)==1
    etaVol=trapz(y,etai);
end

