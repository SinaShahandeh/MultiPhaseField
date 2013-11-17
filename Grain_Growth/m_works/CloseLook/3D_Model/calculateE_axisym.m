% this function calculates energy of the system outside particel in the
% matrix, ME is the Energy density materix with the unit of energy/ volume
% and E is the value of the system energy with the energy unit

function [ME,E]=calculateE_axisym(eta,ppf,mboxsize,nboxsize,delx,settings)

L=settings.L;
alpha=settings.alpha;
beta=settings.beta;
gamma=settings.gamma;
kappa=settings.kappa;
epsilon=settings.epsilon;
DelG=settings.DelG;
accuracy=settings.accuracy;
%% for one-phase field
if size(eta,3)==1
    minE=1/4*alpha*(2*beta-alpha)/beta^2;
    [deletax deletay]=gradient(eta,delx,delx);
    gradeta=sqrt(deletax.^2+deletay.^2);
    % ME=(-alpha/2*eta.^2+beta/4*eta.^4+kappa/2*gradeta.^2+minE);
    ME=(-alpha/2*eta.^2+beta/4*eta.^4 + epsilon*ppf.^2.*eta.^2 + kappa/2*gradeta.^2+minE);

    ME=imcomplement(ppf).*ME;
end

%% for multi-phase field
if size(eta,3)>1
    %% making energy density functional matrix
    minE=1/4*alpha(1)*(2*beta(1)-alpha(1))/beta(1)^2;
    minE=-(-alpha/2+beta/4);
    ME=zeros(mboxsize,nboxsize);
    sumterm=zeros(mboxsize,nboxsize);
    for p=1:size(eta,3)

        try
            [deletax deletay]=gradient(eta(:,:,p),delx,delx);
        catch
            [deletax]=gradient(eta(:,:,p),delx,delx);
            deletay=0;
        end
        gradeta=sqrt(deletax.^2+deletay.^2);
        ME=-alpha/2*eta(:,:,p).^2+beta/4*eta(:,:,p).^4 ...
            +kappa/2*gradeta.^2+ ...
            DelG(p)*eta(:,:,p)+ ...
            ME;
        sumterm=sumterm+(eta(:,:,p).^2.*sum(eta.^2,3)-eta(:,:,p).^4);
    end
    ME=ME+epsilon*ppf.^2.*sum(eta.^2,3)...
        +(gamma*sumterm)+minE;

%removing energy inside the particle
ME=imcomplement(ppf).*ME;
end

%% integration
if strcmp(accuracy,'high')
    x=linspace(0,nboxsize*delx,nboxsize);
    y=linspace(0,mboxsize*delx,mboxsize);
    [X,Y]=meshgrid(x,y);
    E=dblquad(@(xi,yi) interp2(X,Y,ME,xi,yi,'cubic'),0,nboxsize*delx,0,mboxsize*delx,1e-4);
else
    x=linspace(0,nboxsize*delx,nboxsize);
    y=linspace(0,mboxsize*delx,mboxsize);
    [X,Y]=meshgrid(x,y);
    E=(2*pi)^2*sum(sum(ME.*X))*delx^2;
end
