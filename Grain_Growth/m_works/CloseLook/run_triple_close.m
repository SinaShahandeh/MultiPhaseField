%%This is for obtaining equilibrium profiles at triple junctions
clear all
 savedir='/home/cenna/';
% mkdir(savedir)

 figure;
% figure
% phase field parameters
L=5;
alpha=1;
beta=1;
gamma=1;
kappa=2;
epsilon=10;
DelG=[0 0 0];
% setings structure
settings.L=L;
settings.alpha=alpha;
settings.beta=beta;
settings.gamma=gamma;
settings.kappa=kappa;
settings.epsilon=epsilon;
settings.DelG=DelG;
settings.accuracy='low';
% phase fields geometry settings
p=3;
global nboxsize mboxsize
global delx delt scale
scale=2;
mboxsize=23*scale; % y axis in pixels
nboxsize=35*scale; % x axis
delx=2/scale;      % length unit per pixel
grainD=nboxsize+1;  % in pixels

% Particle geometry
pr=3;  % size in length unit
xparticle=fix(nboxsize/2); % position on the grid
yparticle=fix(mboxsize/2*14/10);

endtime=50;
timestepn=30000;
delt=endtime/timestepn;
delt=0.01
% *** Phase Field Procedure *** (so small and simple piece of code!)
%% initial Structure
eta=zeros(mboxsize,nboxsize,p);
% making initial structure
% eta(:,:,1)=circlegrain(mboxsize,nboxsize,nboxsize,6*mboxsize/10,grainD,'dome');
% eta(:,:,2)=circlegrain(mboxsize,nboxsize,0,6*mboxsize/10,grainD,'dome');
% eta(:,:,3)=imcomplement(eta(:,:,1)+eta(:,:,2));

% make 1 central grain 
% % 
% eta(:,:,1)=circlegrain(mboxsize,nboxsize,nboxsize/2,6*mboxsize/10,nboxsize/2,'dome');
% eta(:,:,2)=circlegrain(mboxsize,nboxsize,0,6*mboxsize/10,nboxsize/2,'dome');
% eta(:,:,2)=eta(:,:,2)+circlegrain(mboxsize,nboxsize,nboxsize,6*mboxsize/10,nboxsize/2,'dome');
% eta(:,:,3)=imcomplement(eta(:,:,1)+eta(:,:,2));


% % triple junction equilibriation!
eta(1:mboxsize/3,:,1)=ones(mboxsize/3,nboxsize);
eta(mboxsize/3+1:mboxsize,1:nboxsize/2,2)=ones(mboxsize-mboxsize/3,nboxsize/2);
eta(mboxsize/3+1:mboxsize,nboxsize/2+1:nboxsize,3)=ones(mboxsize-mboxsize/3,nboxsize/2);
% 


% ppf is phase variable representing particles
%particledistro(nboxsize,mboxsize,particles_number,radius)
[ppf]=particledistro(nboxsize,mboxsize,1,pr/delx*2,xparticle,yparticle);

eta2=zeros(mboxsize,nboxsize,p); %pre-assignment

%savesettings
save(strcat(savedir,'/','settings.mat'))

%initialization
for tn=1:5
    for i=1:mboxsize
        for j=1:nboxsize
            del2=1/delx^2*(0.5*(eta(indgsym(i+1,mboxsize),j,:)-2*eta(i,j,:)+eta(indgsym(i-1,mboxsize),j,:))...
                +0.25*(eta(indgsym(i+2,mboxsize),j,:)-2*eta(i,j,:)+eta(indgsym(i-2,mboxsize),j,:)))...
                +1/delx^2*(0.5*(eta(i,indgperi(j+1,nboxsize),:)-2*eta(i,j,:)+eta(i,indgperi(j-1,nboxsize),:))...
                +0.25*(eta(i,indgperi(j+2,nboxsize),:)-2*eta(i,j,:)+eta(i,indgperi(j-2,nboxsize),:)));
            sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
            detadtM=(-alpha.*reshape(eta(i,j,:),1,p)+beta.*reshape(eta(i,j,:),1,p).^3-kappa.*reshape(del2,1,p)+...
                2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j));
            detadt=-L.*(detadtM+2*gamma(1)*reshape(sumterm,1,p));
            eta2(i,j,:)=eta(i,j,:)+reshape(delt*detadt,1,1,p);
            for pind=1:p
                if eta2(i,j,pind)>1
                    eta2(i,j,pind)=1;
                end
                if eta2(i,j,pind)<0
                    eta2(i,j,pind)=0;
                end
            end
        end
    end
    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);
    phi=phi+ppf;
end
se=strel('square',1);
timevec=0;
tn=0;
while tn<timestepn
    Mdetadt=zeros(mboxsize,nboxsize);
    tn=tn+1;
    [yii,xjj]=find(...
        imerode((phi>0.99999),se)==0);
    for ii=1:length(xjj)
        i=yii(ii);j=xjj(ii);
        del2=1/delx^2*(0.5*(eta(indgsym(i+1,mboxsize),j,:)-2*eta(i,j,:)+eta(indgsym(i-1,mboxsize),j,:))...
                +0.25*(eta(indgsym(i+2,mboxsize),j,:)-2*eta(i,j,:)+eta(indgsym(i-2,mboxsize),j,:)))...
                +1/delx^2*(0.5*(eta(i,indgperi(j+1,nboxsize),:)-2*eta(i,j,:)+eta(i,indgperi(j-1,nboxsize),:))...
                +0.25*(eta(i,indgperi(j+2,nboxsize),:)-2*eta(i,j,:)+eta(i,indgperi(j-2,nboxsize),:)));
        sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
        detadtM=(-alpha.*reshape(eta(i,j,:),1,p)+beta.*reshape(eta(i,j,:),1,p).^3-kappa.*reshape(del2,1,p)+...
            2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j)); 
%             +2*ppf(i,j)*[eta(i,j,2)*eta(i,j,3) eta(i,j,1)*eta(i,j,3) eta(i,j,1)*eta(i,j,2)];
        
        detadt=-L.*(detadtM+2*gamma(1)*reshape(sumterm,1,p));
        eta2(i,j,:)=eta(i,j,:)+reshape(delt*detadt,1,1,p);

        for pind=1:p
            if eta2(i,j,pind)>1
                eta2(i,j,pind)=1;
            end
            if eta2(i,j,pind)<0
                eta2(i,j,pind)=0;
            end
        end
    end
    
    % applying boundary condition
%     [ppf2]=particledistro(nboxsize,mboxsize,1,pr/delx*2,xparticle,yparticle);
%     eta2(:,:,1)=boundarycond(eta2(:,:,1),ppf2,0);
%     eta2(:,:,2)=boundarycond(eta2(:,:,2),ppf2,0);

    Mdetadt=(eta2(:,:,1)-eta(:,:,1))/delt;
    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);

    timevec(tn+1)=timevec(tn)+delt;
    %%% Visulaizations

    %% simple structure
    %   subplot(2,3,1)
    %
%     drawgrains(phi,xparticle,yparticle,tn)
%     imshow(eta)
    %     contourgrains(eta,xparticle,yparticle,tn,ppf)
    %% Energy plots
        [ME,E]=calculateE(eta,ppf,mboxsize,nboxsize,delx,settings);
        MME(tn)=E;
    % VecE(tn)=E;
    %     drawE(rot90(rot90(ME)),xparticle,yparticle,tn,eta,ppf)

% 
% 
%     %% Velocity map
%     [nablaetax,nablaetay]=gradient(eta(:,:,1),delx,delx);
%     %     nablaetax((phi>0.99))=nan;
%     %     nablaetay((phi>0.99))=nan;
%     Mvelocx=Mdetadt./nablaetax;
%     Mvelocy=Mdetadt./nablaetay;
%     %         Mveloc=sqrt(Mvelocx.^2+Mvelocy.^2);
%     Mveloc=Mdetadt./sqrt(nablaetax.^2+nablaetay.^2);
%     %
%     %% Average velocity field over certain area of the domain:
%     [vel]=aveveloc(Mveloc,phi,ppf,Mdetadt,'dome');
%     Mvel(tn)=vel;
%     %     MD(tn)=D;
%     h=subplot(2,3,1);
%     plot(timevec(2:end),-Mvel)
%     title('Velocity average based on the field')
%     xlabel(strcat('Time= ', num2str(tn*delt)))
% 
%     %% plotting places where average velocity is calculated
%     maxM=0.9*max(max(abs(Mdetadt)));
%     h=subplot(2,3,6);
%     phivel=phi;
%     phivel(abs(Mdetadt)>maxM)=0;
%     imshow(phivel)
% 
%     %     h=subplot(2,3,2);
%     %     drawvelocity(Mvelocx);
%     %     set(h,'clim',[-0.5 0.5])
%     %     title('X velocity component')
%     %     %
%     %     h=subplot(2,3,3);
%     %     drawvelocity(Mvelocy)
%     %     set(h,'clim',[-0.7 0.7])
%     %     title('Y velocity component')
%     %
%     h=subplot(2,3,2);
%     drawvelocity(Mveloc)
%     set(h,'clim',[-0.9 0])
%     title('Velocity Magnitude')
%     xlabel(strcat('Time= ', num2str(tn*delt)))
%     %
%     %     h=subplot(2,3,6);
%     %     drawvelocity(nablaetax)
%     %     set(h,'clim',[-0.3 0.3])
%     %     title('\nabla \eta x')
%     %
%     %     h=subplot(2,3,4);
%     %     drawvelocity(nablaetay)
%     %     set(h,'clim',[-0.3 0.3])
%     %     title('\nabla \eta y')
%     %
%     %     h=subplot(2,3,1);
%     %     drawvelocity(Mdetadt)
%     %     set(h,'clim',[-0.2 0])
%     %     title('d\eta / dt')
% 
%     %% Force on a circular region inside the domain
%     [F,Fprofile,ME,xp,yp]=Force(eta,ppf,settings,pr,xparticle,yparticle);
%     Mforce(:,tn)=F';
%     h=subplot(2,3,4);
%     plot(timevec(2:end),Mforce(2,:))
%     hold on
%     plot(timevec(2:end),Mforce(1,:),'r')
%     title('Force acting on circular region, red=Fx, blue=Fy')
%     xlabel(strcat('Time= ', num2str(tn*delt)))
%     hold off

    %% Draw Energy Field
%     h=subplot(2,3,3);
    drawvelocity(ME)
    hold on
%     plot3(xp,yp,zeros(1,length(xp))+max(max(ME)),'y')
    % %     set(h,'clim',[-0.2 0])
    title(['Energy density, timestep= ' num2str(tn)])
    hold off
% 
%     %% Volume of each phase field
%     etaVol=etaVolume(eta,delx,nboxsize,mboxsize);
%     MetaVol(:,tn)=etaVol';
%     h=subplot(2,3,5);
%     plot(timevec(2:end),MetaVol(1,:)/(mboxsize*nboxsize*delx^2))
%     hold on
%     plot(timevec(2:end),MetaVol(2,:)/(mboxsize*nboxsize*delx^2),'r')
%     plot(timevec(2:end),(MetaVol(1,:)+MetaVol(2,:))/(mboxsize*nboxsize*delx^2),'g')
%     title('Volume fraction of phases')
%     xlabel(strcat('Time= ', num2str(tn*delt)))
%     hold off
% 


    %   savegrains(eta,ppf,E,xparticle,yparticle,tn,savedir)

    % additional expriments
    %     MMdetadt(:,:,tn)=Mdetadt;
    %     MMeta1(:,:,tn)=eta(:,:,1);
    %
    %% SAVING
% if mod(tn,10)==0
%      save(strcat(savedir,'/',num2str(tn),'.mat'),...
%          'ME','E','eta','tn')
% end
     pause(0.01)
%     tn
    %     disp(strcat('Time= ', num2str(tn*delt)))
end
%

%% 
