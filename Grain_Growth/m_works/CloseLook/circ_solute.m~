%function [InterfaceVel]=realtime_particle_friction_solute(param,filenum,savedir)

clear
for ri=10:10:10
savedir='/home/cenna/Results/2Dclose/circ_solute/mat/'
filenum=0
accuracy='low';
mkdir(savedir)
% phase field parameters
L=1;
m=2;
gamma=1.5*m;
kappa=4;

% solte drag parameters
a=0.0;
b=4;
% settings structural element
settings.L=L(1);
settings.alpha=m;
settings.beta=m;
settings.gamma=gamma(1);
settings.kappa=kappa(1);
settings.epsilon=nan;
settings.DelG=[0 0];

% geometry settings
p=2;
global nboxsize mboxsize
global delx delt scale
scale=2;

delx=2/scale;      % length unit per pixel
r=fix(ri/delx)*scale;
mboxsize=fix((r)+10/delx); % y axis in pixels
nboxsize=fix((r)+10/delx); % x axis

endtime=50;
timestepn=500;
delt=0.05

% *** Phase Field Procedure *** (so small and simple piece of code!)
eta=zeros(mboxsize,nboxsize,p);
% making initial structure
eta(:,:,1)=circlegrain(mboxsize,nboxsize,nboxsize/2,mboxsize/2,r,'circ');
eta(:,:,2)=imcomplement(eta(:,:,1));
eta2=eta;
%savesettings
% save(strcat(pwd,'/',savedir,'/','setings.mat'))

%initialization
for tn=1:50
    tic
    for i=1:mboxsize
        for j=1:nboxsize
            del2=1/delx^2*(0.5*(eta(indg(i+1,mboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-1,mboxsize),j,:)))...
                +1/delx^2*(0.5*(eta(i,indg(j+1,nboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-1,nboxsize),:)));
            sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
            detadtM=m*(-eta(i,j,:)+eta(i,j,:).^3)+2*gamma*sumterm-kappa.*del2;
            detadt=-L.*(detadtM);
            Mdetadt(i,j)=detadt(1);
            eta2(i,j,:)=eta(i,j,:)+delt*detadt;
        end
    end
    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);
%     drawgrains(phi,nan,nan,tn)
    toc
end

se=strel('square',5);
timevec=0;
tn=0;
etaVol=etaVolume(eta(:,:,1),delx,nboxsize,mboxsize,accuracy)
while etaVol>10 %tn<timestepn
    %     Mdetadt=zeros(mboxsize,nboxsize);
    tn=tn+1;
    [yii,xjj]=find(...
        imerode((phi>0.99999),se)==0);

    % Velocity map
    % Mdetadt=(eta2(:,:,1)-eta(:,:,1))/delt;
    [nablaetax,nablaetay]=gradient(eta(:,:,1),delx,delx);
    
    Mveloc=Mdetadt./sqrt(nablaetax.^2+nablaetay.^2);
    % remove wrong data from Mveloc
    [maxdetadt,indimaxdetadt]=max(abs(Mdetadt));
    [maxdetadt,indjmaxdetadt]=max(abs(maxdetadt));
    indimaxdetadt=indimaxdetadt(indjmaxdetadt); %indeceis of maximum in detadt where we know maximum of velocity is calculated correctly
    Mveloc(abs(Mveloc)>abs(Mveloc(indimaxdetadt,indjmaxdetadt)))=nan;
    Mveloc=abs(Mveloc);
    Pf=a*Mveloc./(1+b^2*Mveloc.^2);
    Pf(isnan(Pf))=0;
    for ii=1:length(xjj)
        i=yii(ii);j=xjj(ii);
        del2=1/delx^2*(0.5*(eta(indg(i+1,mboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-1,mboxsize),j,:)))...
                +1/delx^2*(0.5*(eta(i,indg(j+1,nboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-1,nboxsize),:)));
            sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
            detadtM=m*(-eta(i,j,:)+eta(i,j,:).^3)+2*gamma*sumterm-kappa.*del2;
        detadt=-L*(detadtM+6*eta(i,j,:).*(1-eta(i,j,:))*(Pf(i,j))/2);
        % matrix of the eta.dot
        eta2(i,j,:)=eta(i,j,:)+delt*detadt;
        %         Mdel2(i,j)=eta(i,j,2)*del2(1)^2-eta(i,j,1)*del2(2)^2;
        Mdetadt(i,j)=detadt(1);
    end
    %     Mdetadt=(eta2(:,:,1)-eta(:,:,1))/delt;
    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);
    timevec(tn+1)=timevec(tn)+delt;
    %%% Visulaizations

    %% simple structure
    %   subplot(2,3,1)
    %
%      drawgrains(phi,nan,nan,tn)
    %     contourgrains(eta,xparticle,yparticle,tn,ppf)
    %% Energy plots
    %     [ME,E]=calculateE(eta,ppf,mboxsize,nboxsize,delx);
    % VecE(tn)=E;
    %     drawE(rot90(rot90(ME)),xparticle,yparticle,tn,eta,ppf)

    %% Middle point curvature and speed
    %     [cur1,cur2,detach,cent1,cent2]=analysecontourgrains(eta,xparticle,yparticle,tn,ppf);
    %     curi1(tn)=cur1;
    %     curi2(tn)=cur2;
    %     detachi(tn)=detach;
    %% speed
    %     speed1(tn)=(cent1-cent1p)/delt;
    %     speed2(tn)=(cent2-cent2p)/delt;
    %     cent1p=cent1;
    %     cent2p=cent2;
    %     Mcent1(tn)=cent1;Mcent2(tn)=cent2;

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



    %% Average velocity field over certain area of the domain:
    %         [vel]=aveveloc(Mveloc,phi,ppf,Mdetadt,'dome');
    %         Mvel(tn)=vel;
    %         %     MD(tn)=D;
    %         h=subplot(2,2,1);
    %         plot(timevec(2:end),-Mvel)
    %         title('Velocity average based on the field')
    %         xlabel(strcat('Time= ', num2str(tn*delt)))

    %% plotting places where average velocity is calculated
    %     maxM=0.9*max(max(abs(Mdetadt)));
    %     h=subplot(2,3,6);
    %     phivel=phi;
    %     phivel(abs(Mdetadt)>maxM)=0;
    %     imshow(phivel)

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
    %     h=subplot(2,3,2);
    %     drawvelocity(Mveloc)
    %     set(h,'clim',[-0.9 0])
    %     title('Velocity Magnitude')
    %     xlabel(strcat('Time= ', num2str(tn*delt)))
    %     %
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
    %     [F,Fprofile,ME,xp,yp]=Force(eta,ppf,setings,pr,xparticle,yparticle);
    %     Mforce(:,tn)=F';
    %     h=subplot(2,3,4);
    %     plot(timevec(2:end),Mforce(2,:))
    %     hold on
    %     plot(timevec(2:end),Mforce(1,:),'r')
    %     title('Force acting on circular region, red=Fx, blue=Fy')
    %     xlabel(strcat('Time= ', num2str(tn*delt)))
    %     hold off

    %% Draw Energy Field
    %         h=subplot(2,2,2);
    %     [ME,E]=calculateE(eta,ppf,mboxsize,nboxsize,delx,settings);
    %     drawvelocity(ME)
    % InterfaceVel
    %     % %     set(h,'clim',[-0.2 0])
    %     title(['Energy density at timestep' num2str(tn)])


    %% Volume of each phase field
    etaVol=etaVolume(eta(:,:,1),delx,nboxsize,mboxsize,accuracy);
    MetaVol(:,tn)=etaVol;
    %         h=subplot(2,2,3);
    %         plot(timevec(2:end),MetaVol(1,:)/(mboxsize*nboxsize*delx^2))
    %         hold on
    %         plot(timevec(2:end),MetaVol(2,:)/(mboxsize*nboxsize*delx^2),'r')
    %         plot(timevec(2:end),(MetaVol(1,:)+MetaVol(2,:))/(mboxsize*nboxsize*delx^2),'g')
    %         title('Volume fraction of phases')
    %         xlabel(strcat('Time= ', num2str(tn*delt)))
    %         hold off

    %% Draw del2 Field
    %    h=subplot(2,1,1);
    %       drawvelocity(Mdel2)
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

    %     disp(strcat('Time= ', num2str(tn*delt)))
end
%


save([savedir num2str(filenum) '.mat'])

%% Velocity for circle
% 
mobility=3/2*L(1)*sqrt(2*kappa(1)/m);
intenergy=1/3*sqrt(2*m*kappa);
 time=timevec;time(1)=[];
    vol=MetaVol;
    radius=sqrt(vol/pi);
    n=1;steps=[1];
%     plot(time,radius)
    spline1 = spaps(time(1:steps(n):end),radius(1:steps(n):end),0.000001,3);
    Mvel=fnval(fnder(spline1,1),time(1:steps(n):end));
% xlabel('Time')
% ylabel('radius of circular grain')
% figure
Mvel=-Mvel;
Mcurvature=1./radius(1:steps(n):end);
plot(Mcurvature,Mvel/mobility/intenergy)
xlabel('\Delta G/ \sigma_{gb}');
ylabel('V / (M \sigma_{gb})');
grid on

axis([0 0.5 0 0.5]);box on
hold on
plot([0 0.5],[0 0.5],'r')


%% calculating average speed for dome
% gamma=1.5;
% fsaddle=(2*gamma-1)/(4*(2*gamma+1));
% intenergy=8/3*fsaddle*sqrt(2*m*kappa);
% mobility=3/2*L*sqrt(2*kappa/m);
% % mobility=3/2*L*sqrt(2*kappa/m);
% % intenergy=1/3*sqrt(2*m*kappa);
% w=quad(@(x) interp1(delx*[1:nboxsize], eta(mboxsize,:,1),x,'spline'),0,delx*nboxsize);
% pp=polyfit(timevec(1:end-1),-MetaVol(1,:)/w,1);
% InterfaceVel=pp(1)
% Intcurvature=1./(r/2*delx);
% hold on
% plot(Intcurvature, InterfaceVel/mobility/intenergy,'o')
% xlabel('\Delta G/ \sigma_{gb}');
% ylabel('V / (M \sigma_{gb})');
% grid on
% hold on
% plot([0 0.5],[0 0.5],'r')
end