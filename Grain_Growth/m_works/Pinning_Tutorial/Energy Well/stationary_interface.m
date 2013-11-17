%% this is for straigth line diffuse interfacer with 2 phase field. and a


% function [eta,ppf,ME,E]=stationary_interface(x,epsilon,pr)
% Initial position of interface from bottom of the domain (x=from 0 to 1)
x=0.5;
% Particle radius in unit of x
pr=8;

% figure;
% phase field parameters
L=5*[1 1];
alpha=[1 1];
beta=[1 1];
gamma=1.5*[1 1];
kappa=2*[1 1];
epsilon=5;
DelG=[0 0];
% BCValue=value;

% setings structure
settings.L=L(1);
settings.alpha=alpha(1);
settings.beta=beta(1);
settings.gamma=gamma(1);
settings.kappa=kappa(1);
settings.epsilon=epsilon(1);
settings.DelG=DelG;
settings.accuracy='low';
% geometry settings
p=2;
global nboxsize mboxsize
global delx delt scale

% Domain Size
scale=3;           % larger scale makes delx smaller without changing the system size
mboxsize=35*scale; % y axis in pixels
nboxsize=35*scale; % x axis
delx=2/scale;     

% time settings
timestepn=100;
delt=0.01;

% Particle geometry
% pr=3;  % size in length unit
xparticle=fix(nboxsize/2); % position on the grid (pixel)
yparticle=fix(mboxsize/2*10/10);
% ppf is phase variable representing particles
%particledistro(nboxsize,mboxsize,particles_number,radius)
[ppf]=particledistro(nboxsize,mboxsize,1,pr/delx*2,xparticle,yparticle);
eta2=zeros(mboxsize,nboxsize,p); %pre-assignment

% *** Phase Field Procedure *** (so small and simple piece of code!)
eta=zeros(mboxsize,nboxsize,p);
% making initial structure

eta(:,:,1)=zeros(mboxsize,nboxsize);
x=fix(mboxsize*x);
eta(1:x,:,1)=1;
eta(:,:,2)=imcomplement(eta(:,:,1));

%initialization
for tn=1:5
    for i=1:mboxsize
        for j=1:nboxsize
            del2=1/delx^2*(0.5*(eta(indg(i+1,mboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-1,mboxsize),j,:))...
                +0.25*(eta(indg(i+2,mboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-2,mboxsize),j,:)))...
                +1/delx^2*(0.5*(eta(i,indg(j+1,nboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-1,nboxsize),:))...
                +0.25*(eta(i,indg(j+2,nboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-2,nboxsize),:)));
            sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
            detadtM=(-alpha.*reshape(eta(i,j,:),1,p)+beta.*reshape(eta(i,j,:),1,p).^3-kappa.*reshape(del2,1,p)+...
                2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j));
            detadt=-L.*(detadtM+2*gamma(1)*reshape(sumterm,1,p));
            eta2(i,j,:)=eta(i,j,:)+reshape(delt*detadt,1,1,2);
            for pind=1:p
                if eta2(i,j,pind)>1
                    eta2(i,j,pind)=1;
                end
                if eta2(i,j,pind)<0
                    eta2(i,j,pind)=0;
                end
            end
            Mdel2(i,j)=del2(2);
        end
    end
    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);
    phi=phi+ppf;
end
se=strel('square',6);
timevec=0;
tn=0;
while tn<timestepn
    Mdetadt=zeros(mboxsize,nboxsize);
    tn=tn+1;
    [yii,xjj]=find(...
        imerode((phi>0.9),se)==0);
    for ii=1:length(xjj)
        i=yii(ii);j=xjj(ii);
        del2=1/delx^2*(0.5*(eta(indg(i+1,mboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-1,mboxsize),j,:))...
            +0.25*(eta(indg(i+2,mboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-2,mboxsize),j,:)))...
            +1/delx^2*(0.5*(eta(i,indg(j+1,nboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-1,nboxsize),:))...
            +0.25*(eta(i,indg(j+2,nboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-2,nboxsize),:)));
        sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
        detadtM=(-alpha.*reshape(eta(i,j,:),1,p)+beta.*reshape(eta(i,j,:),1,p).^3-kappa.*reshape(del2,1,p)+...
            2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j));
        detadt=-L.*(detadtM+2*gamma(1)*reshape(sumterm,1,p));

        eta2(i,j,:)=eta(i,j,:)+reshape(delt*detadt,1,1,2);

        for pind=1:p
            if eta2(i,j,pind)>1
                eta2(i,j,pind)=1;
            end
            if eta2(i,j,pind)<0
                eta2(i,j,pind)=0;
            end
        end
        Mdel2(i,j)=del2(2);
    end


    eta=eta2;
    % making optimization matrix
    [deletax1 deletay1]=gradient(eta(:,:,1),delx,delx);
    [deletax2 deletay2]=gradient(eta(:,:,2),delx,delx);
    phi=imcomplement(abs(deletax1)+abs(deletay1)+abs(deletax2)+abs(deletay2));


    timevec(tn+1)=timevec(tn)+delt;
    %%% Visulaizations

    %% Draw Energy Field
    [ME,E]=calculateE(eta,ppf,mboxsize,nboxsize,delx,settings);
    drawvelocity(ME)
    title(['Energy density at timestep' num2str(tn)])
    % %     set(h,'clim',[-0.2 0])
    pause(0.01)
end

disp('Total energy of the sytem:')
E

