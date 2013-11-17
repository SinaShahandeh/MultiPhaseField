% this function calculates energy of the system outside particel in the
% matrix, ME is the Energy density materix with the unit of energy/ volume
% and E is the value of the system energy with the energy unit

function [ME,E]=calculateE_3D(eta,ppf,delx,settings)

L=settings.L;
alpha=settings.alpha;
beta=settings.beta;
gamma=settings.gamma;
kappa=settings.kappa;
epsilon=settings.epsilon;
DelG=settings.DelG;
accuracy=settings.accuracy;

mboxsize=size(eta,1);
nboxsize=size(eta,2);
lboxsize=size(eta,3);
%% for multi-phase field
if size(eta,3)>1
    %% making energy density functional matrix
    minE=1/4*alpha(1)*(2*beta(1)-alpha(1))/beta(1)^2;
    minE=-(-alpha/2+beta/4);
    ME=zeros(mboxsize,nboxsize,lboxsize);
    sumterm=zeros(mboxsize,nboxsize,lboxsize);
    for p=1:size(eta,4)
        [deletax deletay deletaz]=gradient(eta(:,:,:,p),delx,delx,delx);
        gradeta=sqrt(deletax.^2+deletay.^2+deletaz.^2);
        ME=-alpha/2*eta(:,:,:,p).^2+beta/4*eta(:,:,:,p).^4 ...
            +kappa/2*gradeta.^2+ ...
            DelG(p)*eta(:,:,:,p)+ ...
            ME;
        sumterm=sumterm+(eta(:,:,:,p).^2.*sum(eta.^2,4)-eta(:,:,:,p).^4);
    end
    ME=ME+epsilon*ppf.^2.*sum(eta.^2,4)...
        +(gamma*sumterm)+minE;
    % removing energy inside the particle
    ME=imcomplement(ppf).*ME;
end

%% integration
if strcmp(accuracy,'high')
    x=linspace(0,nboxsize*delx,nboxsize);
    y=linspace(0,mboxsize*delx,mboxsize);
    [X,Y]=meshgrid(x,y);
    E=dblquad(@(xi,yi) interp2(X,Y,ME,xi,yi,'cubic'),0,nboxsize*delx,0,mboxsize*delx,1e-4);
else
    E=sum(sum(sum(ME)))*delx^3;
end
