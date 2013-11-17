%% Force acting on the specific path on the multi phase field domain

function [F,Fprofile,ME,xp,yp]=Force(eta,ppf,setings,pr,xparticle,yparticle)
global delx scale
global mboxsize nboxsize

alpha=setings.alpha;
beta=setings.beta;
gamma=setings.gamma;
kappa=setings.kappa;
epsilon=setings.epsilon;

%% making energy density functional matrix
minE=1/4*alpha(1)*(2*beta(1)-alpha(1))/beta(1)^2;
ME=zeros(mboxsize,nboxsize);
for p=1:size(eta,3)
    
    [deletax deletay]=gradient(eta(:,:,p),delx,delx);
    gradeta=sqrt(deletax.^2+deletay.^2);
    ME=-alpha(p)/2*eta(:,:,p).^2+beta(p)/4*eta(:,:,p).^4 ...
    +kappa(p)/2*gradeta.^2+ ...
    (gamma(p)*eta(:,:,p).^2.*sum(eta,3).^2-eta(:,:,p).^4) ...
    +minE+ME;
end
ME=ME+epsilon*ppf.^2.*sum(eta.^2,3);

% removing energy inside the particle
% ME=imcomplement(ppf).*ME;

%% performing line integral of the energy field


diameter=2*pr/delx; %in grid point

t=linspace(0,2*pi,1000);


if mod(diameter,2)==0
    y1=yparticle-fix(diameter/2)+1;y2=yparticle+fix(diameter/2);
    x1=xparticle-fix(diameter/2)+1;x2=xparticle+fix(diameter/2);
else
    y1=yparticle-fix(diameter/2);y2=yparticle+fix(diameter/2);
    x1=xparticle-fix(diameter/2);x2=xparticle+fix(diameter/2);
end
if y1<1
    y1=1;
end
if y2>mboxsize
    y2=mboxsize;
end
if x1<1
    x1=1;
end
if x2>nboxsize
    x2=nboxsize;
end


% defining the path
% go a little bit furthur from the particle radius
pr=pr+1*pr;
xp=xparticle*delx+pr*cos(t);
yp=yparticle*delx+pr*sin(t);
% integral of the path length s=int(ds)=int(pr*dt)
s=pr*t;

%% evaluating energy density along the path 
xi=[0:nboxsize-1]*delx;
yi=[0:mboxsize-1]*delx;
[X,Y]=meshgrid(xi,yi);
f=interp2(X,Y,ME,xp,yp,'spline');

%the normal vwctor angle in here is equal to the t (fortunately!)
fx=f.*cos(t);
fy=f.*sin(t);

Fprofile=[xp;yp;fx;fy];

F=[trapz(s,fx) trapz(s,fy)];