clear all
savedir='velocity1';
mkdir(savedir)
% figure;
clf
% phase field parameters
L=5*[1 1];
alpha=[1 1];
beta=[1 1];
gamma=[1 1];
kappa=[1 1];
epsilon=5;
G=[0 0];
% setings structure
setings.L=L;
setings.alpha=alpha;
setings.beta=beta;
setings.gamma=gamma;
setings.kappa=kappa;
setings.epsilon=epsilon;

% geometry settings
p=2;
global nboxsize mboxsize
global delx delt scale
scale=3;
mboxsize=23*scale; % y axis in pixels
nboxsize=22*scale; % x axis
delx=2/scale;      % length unit per pixel
grainD=15*scale;  % in pixels

% Particle geometry
pr=1;  % size in length unit
xparticle=fix(nboxsize/2); % position on the grid
yparticle=fix(mboxsize/2*7.5/10);

endtime=50;
timestepn=3000;
delt=endtime/timestepn;


% *** Phase Field Procedure *** (so small and simple piece of code!)
eta=zeros(mboxsize,nboxsize,p);
% making initial structure
eta(:,:,1)=circlegrain(mboxsize,nboxsize,nboxsize/2,4.5*mboxsize/10,grainD,'dome');
eta(:,:,2)=imcomplement(eta(:,:,1));
% ppf is phase variable representing particles
%particledistro(nboxsize,mboxsize,particles_number,radius)
[ppf]=particledistro(nboxsize,mboxsize,0,pr/delx*2,xparticle,yparticle);

eta2=zeros(mboxsize,nboxsize,p); %pre-assignment

%savesettings
save(strcat(pwd,'/',savedir,'/','setings.mat'))

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
                2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j))+G;
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
        end
    end
    
    [ppf2]=particledistro(nboxsize,mboxsize,1,pr/delx*2,xparticle,yparticle);
    eta2(:,:,1)=boundarycond(eta2(:,:,1),ppf2,0);
    eta2(:,:,2)=boundarycond(eta2(:,:,2),ppf2,0);
   
    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);
    phi=phi+ppf;
    [cur1,cur2,detach,cent1,cent2]=analysecontourgrains(eta,xparticle,yparticle,tn,ppf);
    cent1p=cent1;
    cent2p=cent2;
end
se=strel('square',5);
timevec=0;
tn=0;
while tn<3000
    Mdetadt=zeros(mboxsize,nboxsize);
    tn=tn+1;
    [yii,xjj]=find(...
        imerode((phi>0.99999),se)==0);
    for ii=1:length(xjj)
        i=yii(ii);j=xjj(ii);
        del2=1/delx^2*(0.5*(eta(indg(i+1,mboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-1,mboxsize),j,:))...
            +0.25*(eta(indg(i+2,mboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-2,mboxsize),j,:)))...
            +1/delx^2*(0.5*(eta(i,indg(j+1,nboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-1,nboxsize),:))...
            +0.25*(eta(i,indg(j+2,nboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-2,nboxsize),:)));
        sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
        detadtM=(-alpha.*reshape(eta(i,j,:),1,p)+beta.*reshape(eta(i,j,:),1,p).^3-kappa.*reshape(del2,1,p)+...
            2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j))+G;
        detadt=-L.*(detadtM+2*gamma(1)*reshape(sumterm,1,p));
        % matrix of the eta.dot
        eta2(i,j,:)=eta(i,j,:)+reshape(delt*detadt,1,1,2);
        %             elovution of particle
        %             del2=1/delx^2*(0.5*(ppf(indg(i+1,gridn),j)-2*ppf(i,j)+ppf(indg(i-1,gridn),j))...
        %                 +0.25*(ppf(indg(i+2,gridn),j)-2*ppf(i,j)+ppf(indg(i-2,gridn),j)))...
        %                 +1/delx^2*(0.5*(ppf(i,indg(j+1,gridn))-2*ppf(i,j)+ppf(i,indg(j-1,gridn)))...
        %                 +0.25*(ppf(i,indg(j+2,gridn))-2*ppf(i,j)+ppf(i,indg(j-2,gridn))));
        %             sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
        %             detadtM=(-alpha.*reshape(eta(i,j,:),1,p)+beta.*reshape(eta(i,j,:),1,p).^3-kappa.*reshape(del2,1,p)+...
        %                 2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j))+G;
        %             detadt=-L.*(detadtM+2*gamma*reshape(sumterm,1,p));
        %             eta2(i,j,:)=eta(i,j,:)+reshape(delt*detadt,1,1,2);
        for pind=1:p
            if eta2(i,j,pind)>1
                eta2(i,j,pind)=1;
            end
            if eta2(i,j,pind)<0
                eta2(i,j,pind)=0;
            end
        end
    end
    Mdetadt=(eta2(:,:,1)-eta(:,:,1))/delt;
    
    
    [ppf2]=particledistro(nboxsize,mboxsize,1,pr/delx*2,xparticle,yparticle);
    eta2(:,:,1)=boundarycond(eta2(:,:,1),ppf2,0);
    eta2(:,:,2)=boundarycond(eta2(:,:,2),ppf2,0);
    
    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);
    
    timevec(tn+1)=timevec(tn)+delt;
    %%% Visulaizations

    %% simple structure
    %   subplot(2,3,1)
    %
    % drawgrains(phi,xparticle,yparticle,tn)
    %     contourgrains(eta,xparticle,yparticle,tn,ppf2)
    %% Energy plots
    %     [ME,E]=calculateE(eta,ppf2,mboxsize,nboxsize,delx);
    % VecE(tn)=E;
    %     drawE(rot90(rot90(ME)),xparticle,yparticle,tn,eta,ppf2)

    %% Middle point curvature and speed
    [cur1,cur2,detach,cent1,cent2]=analysecontourgrains(eta,xparticle,yparticle,tn,ppf2);
    curi1(tn)=cur1;
    curi2(tn)=cur2;
    detachi(tn)=detach;
    %% speed
    %     speed1(tn)=(cent1-cent1p)/delt;
    %     speed2(tn)=(cent2-cent2p)/delt;
    %     cent1p=cent1;
    %     cent2p=cent2;
    Mcent1(tn)=cent1;Mcent2(tn)=cent2;

    %
    %     %% Speed Plots
    %     subplot(2,3,[6])
    %     plot(timevec(2:end),speed1,timevec(2:end),speed2,'r');
    %     title(strcat('v_{\eta_1}= ', num2str(speed1(tn)),' , ',...
    %         'v_{\eta_2}= ', num2str(speed2(tn))));
    %     xlabel('time')
    %     %% Curvature Plots
    %     subplot(3,1,3)
    %     plot(timevec(2:end),curi1,timevec(2:end),curi2,'r');
    %         xlabel('time')
    %     title(strcat('\kappa_{\eta_1}= ', num2str(cur1),' , ',...
    %         '\kappa_{\eta_2}= ', num2str(cur2)));


    %% Velocity map
    [nablaetax,nablaetay]=gradient(eta(:,:,1),delx,delx);
    %     nablaetax((phi>0.99))=nan;
    %     nablaetay((phi>0.99))=nan;
    Mvelocx=Mdetadt./nablaetax;
    Mvelocy=Mdetadt./nablaetay;
    %         Mveloc=sqrt(Mvelocx.^2+Mvelocy.^2);
    Mveloc=Mdetadt./sqrt(nablaetax.^2+nablaetay.^2);
    %
    %% Average velocity field over certain area of the domain:
    [vel]=aveveloc(Mveloc,phi,ppf2,Mdetadt,'dome');
    Mvel(tn)=vel;
    %     MD(tn)=D;
    h=subplot(2,3,1);
    plot(timevec(2:end),-Mvel)
    title('Velocity average based on the field')
    xlabel(strcat('Time= ', num2str(tn*delt)))

    %% plotting places where average velocity is calculated
    maxM=0.9*max(max(abs(Mdetadt)));
    h=subplot(2,3,6);
    phivel=phi;
    phivel(abs(Mdetadt)>maxM)=0;
    imshow(phivel)

    %     h=subplot(2,3,2);
    %     drawvelocity(Mvelocx);
    %     set(h,'clim',[-0.5 0.5])
    %     title('X velocity component')
    %     %
    %     h=subplot(2,3,3);
    %     drawvelocity(Mvelocy)
    %     set(h,'clim',[-0.7 0.7])
    %     title('Y velocity component')
    %
    h=subplot(2,3,2);
    drawvelocity(Mveloc)
    set(h,'clim',[-0.9 0])
    title('Velocity Magnitude')
    xlabel(strcat('Time= ', num2str(tn*delt)))
    %
    %     h=subplot(2,3,6);
    %     drawvelocity(nablaetax)
    %     set(h,'clim',[-0.3 0.3])
    %     title('\nabla \eta x')
    %
    %     h=subplot(2,3,4);
    %     drawvelocity(nablaetay)
    %     set(h,'clim',[-0.3 0.3])
    %     title('\nabla \eta y')
    %
    %     h=subplot(2,3,1);
    %     drawvelocity(Mdetadt)
    %     set(h,'clim',[-0.2 0])
    %     title('d\eta / dt')

    %% Force on a circular region inside the domain
    [F,Fprofile,ME,xp,yp]=Force(eta,ppf2,setings,pr,xparticle,yparticle);
    Mforce(:,tn)=F';
    h=subplot(2,3,4);
    plot(timevec(2:end),Mforce(2,:))
    hold on
    plot(timevec(2:end),Mforce(1,:),'r')
    title('Force acting on circular region, red=Fx, blue=Fy')
    xlabel(strcat('Time= ', num2str(tn*delt)))
    hold off

    %% Draw Energy Field
    h=subplot(2,3,3);
    drawvelocity(ME)
    hold on
    plot3(xp,yp,zeros(1,length(xp))+max(max(ME)),'y')
    % %     set(h,'clim',[-0.2 0])
    title('Energy density')
    hold off

    %% Volume of each phase field
    etaVol=etaVolume(eta);
    MetaVol(:,tn)=etaVol';
    h=subplot(2,3,5);
    plot(timevec(2:end),MetaVol(1,:)/(mboxsize*nboxsize*delx^2))
    hold on
    plot(timevec(2:end),MetaVol(2,:)/(mboxsize*nboxsize*delx^2),'r')
    plot(timevec(2:end),(MetaVol(1,:)+MetaVol(2,:))/(mboxsize*nboxsize*delx^2),'g')
    title('Volume fraction of phases')
    xlabel(strcat('Time= ', num2str(tn*delt)))
    hold off



    %   savegrains(eta,ppf2,E,xparticle,yparticle,tn,savedir)

    % additional expriments
    %     MMdetadt(:,:,tn)=Mdetadt;
    %     MMeta1(:,:,tn)=eta(:,:,1);
    %
    %% SAVING

    %  save(strcat(pwd,'/',savedir,'/',num2str(tn),'.mat'),...
    %      'Mdetadt','eta','Mveloc','Mvelocx','Mvelocy','tn')

    pause(0.05)
    %     disp(strcat('Time= ', num2str(tn*delt)))
end
%

