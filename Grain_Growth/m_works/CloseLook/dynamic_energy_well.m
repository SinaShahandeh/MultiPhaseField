% function dynamic_energy_well(param)
clear all

x=0.25; % initial condition.
xend=0.85; % ending potions of interface
value=5e-18;
pr=6
particles_number=1
% figure;
% phase field parameters
L=5*[1 1]*1e18;
alpha=1*[1 1]*1e-18;
beta=1*[1 1]*1e-18;
gamma=1*1.5*[1 1]*1e-18;
kappa=0.5*[1 1]*1e-18;
epsilon=value;

DelG=[0 1]*1e-19;
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
scale=2;
mboxsize=30*scale; % y axis in pixels
nboxsize=20*scale; % x axis
delx=2/scale;      % length unit per pixel

% Particle geometry
% pr=3;  % size in length unit
xparticle=fix(nboxsize/2); % position on the grid (pixel)
yparticle=fix(mboxsize/2*10/10);

timestepn=2000;
delt=0.02;
%savesettings
% save(strcat(pwd,'/',savedir,'/','setings.mat'))

% *** Phase Field Procedure *** (so small and simple piece of code!)
eta=zeros(mboxsize,nboxsize,p);
% making initial structure
% eta(:,:,1)=circlegrain(mboxsize,nboxsize,nboxsize/2,4.5*mboxsize/10,grainD,'dome');
eta(:,:,1)=zeros(mboxsize,nboxsize);
x=fix(mboxsize*x);
eta(1:x,:,1)=1;
eta(:,:,2)=imcomplement(eta(:,:,1));
% ppf is phase variable representing particles
%particledistro(nboxsize,mboxsize,particles_number,radius)
[ppf]=particledistro(nboxsize,mboxsize,particles_number,pr/delx*2,xparticle,yparticle);
ppf2=ppf;
eta2=zeros(mboxsize,nboxsize,p); %pre-assignment
%initialization

%initialization
for tn=1:50
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
        end
    end
    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);
    drawgrains(phi,xparticle,yparticle,tn)
    
    %     [cur1,cur2,detach,cent1,cent2]=analysecontourgrains(eta,xparticle,yparticle,tn,ppf);
    %     cent1p=cent1;
    %     cent2p=cent2;
    pause(0.01)
end
se=strel('square',6);
timevec=0;
tn=0;
xend=fix(mboxsize*xend);
while eta(xend,1,1)<0.5 % when interface arrives at xend
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
            2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j))...
            +6*(reshape(eta(i,j,:),1,p)-reshape(eta(i,j,:),1,p).^2).*DelG;
        detadt=-L.*(detadtM+2*gamma(1)*reshape(sumterm,1,p));
        % matrix of the eta.dot
        eta2(i,j,:)=eta(i,j,:)+reshape(delt*detadt,1,1,2);

        %         Mdel2(i,j)=del2(2);
    end

    %     [ppf2]=particledistro(nboxsize,mboxsize,1,pr/delx*2,xparticle,yparticle);
    %     eta2(:,:,1)=boundarycond(eta2(:,:,1),ppf2,0);
    %     eta2(:,:,2)=boundarycond(eta2(:,:,2),ppf2,0);


    %     Mdetadt=(eta2(:,:,1)-eta(:,:,1))/delt;
    eta=eta2;

    % making optimization matrix
    [deletax1 deletay1]=gradient(eta(:,:,1),delx,delx);
    [deletax2 deletay2]=gradient(eta(:,:,2),delx,delx);
    phi=imcomplement(abs(deletax1)+abs(deletay1)+abs(deletax2)+abs(deletay2));
    %     phi=sum(eta(:,:,1:p).^2,3);

    timevec(tn+1)=timevec(tn)+delt;
    %%% Visulaizations

    %% simple structure
    %   subplot(2,3,1)
    %
     drawgrains(phi,xparticle,yparticle,tn)
    %     contourgrains(eta,xparticle,yparticle,tn,ppf)

    %% Draw Energy Field
    %    h=subplot(2,1,2);
    [ME,E]=calculateE(eta,ppf,mboxsize,nboxsize,delx,settings);
     drawvelocity(ME)
    MME(tn)=E;
      title(['Energy density at timestep' num2str(tn)])
    % %     set(h,'clim',[-0.2 0])

    %% Volume of each phase field
    etaVol=etaVolume(eta(:,:,1),delx,nboxsize,mboxsize,'low');
    MetaVol(tn)=etaVol;
    %     h=subplot(2,3,5);
    %     plot(timevec(2:end),MetaVol(1,:)/(mboxsize*nboxsize*delx^2))
    %     hold on
    %     plot(timevec(2:end),MetaVol(2,:)/(mboxsize*nboxsize*delx^2),'r')
    %     plot(timevec(2:end),(MetaVol(1,:)+MetaVol(2,:))/(mboxsize*nboxsize*delx^2),'g')
    %     title('Volume fraction of phases')
    %     xlabel(strcat('Time= ', num2str(tn*delt)))
    %     hold off

    % %% Draw del2 Field
    %    h=subplot(2,1,1);
    %
    %     drawvelocity(Mdel2)
    %
    %     % %     set(h,'clim',[-0.2 0])
    %     title(['\nabla^2 at timestep' num2str(tn)])

    %   savegrains(eta,ppf,E,xparticle,yparticle,tn,savedir)

    % additional expriments
    %     MMdetadt(:,:,tn)=Mdetadt;
    %     MMeta1(:,:,tn)=eta(:,:,1);
    %
    %% SAVING

    %  save(strcat(pwd,'/',savedir,'/',num2str(tn),'.mat'),...
    %      'Mdetadt','eta','Mveloc','Mvelocx','Mvelocy','tn')
    pause(0.01)
    %     disp(strcat('Time= ', num2str(tn*delt)))

    %% find position of interface on the corner
    intindex=find(eta(:,end,1)>0.01 & eta(:,end,1)<0.99);
    intpos=interp1(eta(intindex,end,1),intindex*delx,0.5,'spline');
    Mintpos(tn)=intpos;
        if mod(tn,200)==0
                drawvelocity(ME)
                title(['Energy density at timestep' num2str(tn)])
%               filename=strcat(pwd,'/',num2str(tn),'.png');
%             print('-f1','-dpng','-r100',filename)
        end
end

% % plotting profiles of energy order parameters
% figure
% y=[1:mboxsize]*delx-mboxsize*delx/2;
% subplot(2,1,1)
% plot(y,eta(:,fix(nboxsize/2),1))
% hold on
% plot(y,eta(:,fix(nboxsize/2),2),'r')
% ylabel('Ordere Parameter')
% plot(y,ppf(:,fix(nboxsize/2)),'g')
%
% plot(mboxsize/2*delx*[0 0],[0 1],'k')
%
% subplot(2,1,2)
% plot(y,ME(:,fix(nboxsize/2)),'r')
% xlabel('Position (X unit)')
% ylabel('Energy Density')
% hold on
% plot(mboxsize/2*delx*[0 0],[0 0.8e-18],'k')

% Sigma=E/(nboxsize*delx)*1e18


%% calculating potential well
% figure
% plot(timevec(2:end),MME)

% figure
% plot(Mintpos,MME)
% hold on
% %fitting two straight lines on energy curve
% ln1ind=find(Mintpos>0 & Mintpos<21);
% plot(Mintpos(ln1ind),MME(ln1ind),'y')
% pp1=polyfit(Mintpos(ln1ind),MME(ln1ind),1);
% 
% ln2ind=find(Mintpos>36);
% plot(Mintpos(ln2ind),MME(ln2ind),'g')
% pp2=polyfit(Mintpos(ln2ind),MME(ln2ind),1);
% % potential well re-construction
% well=MME-polyval(pp1,Mintpos);
% figure
% plot(Mintpos,well)

%% Plotting based on equvalent distance comes from volume of the order parameter

figure
width=nboxsize*delx;
Mintpos=MetaVol/width;
plot(MetaVol/width,MME)
hold on
%fitting two straight lines on energy curve
ln1ind=find(Mintpos>15 & Mintpos<16);
plot(MetaVol(ln1ind)/width,MME(ln1ind),'y')
pp1=polyfit(MetaVol(ln1ind)/width,MME(ln1ind),1);

ln2ind=find(Mintpos>40);
plot(MetaVol(ln2ind)/width,MME(ln2ind),'g')
pp2=polyfit(MetaVol(ln2ind)/width,MME(ln2ind),1);

% potential well re-construction

well=MME-polyval((1*pp1+0.0*pp2),MetaVol/width);
% well=MME-polyval(pp1,MetaVol/width);

figure
plot(MetaVol/width,well)
xlabel('Equivalent Interface Position (nm)')
ylabel('Total System Energy (J)')




