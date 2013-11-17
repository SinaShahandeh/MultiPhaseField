clear
savedir='/home/cenna/Results/2D/Fric800s2_m1_k4/0/';
start=2000;
step=2000;
ending=600000;
tn=10000
mboxsize=1600;
nboxsize=mboxsize;
tni=0;
% for tn=[start:step:ending]
    %     tni=tni+1;
    Kdata=importdata([savedir 'K_' num2str(tn) '.txt']);
    Kraw=Kdata(:,4);Kraw(abs(Kraw)>0.04)=[];
    Khist=abs(Kraw);
    
    Kneg=zeros(mboxsize,nboxsize)*nan;
%     Kindsneg=zeros(mboxsize,nboxsize);
%     Kindspos=zeros(mboxsize,nboxsize);
    for i=1:length(Kdata)
        if Kdata(i,4)<=0 % select the negative component
            Kneg(Kdata(i,1)+1,Kdata(i,2)+1)=Kdata(i,4);
%             Kindsneg(Kdata(i,1)+1,Kdata(i,2)+1)=Kdata(i,3);
        end
    end
    Kpos=zeros(mboxsize,nboxsize)*nan;
    for i=1:length(Kdata)
        if Kdata(i,4)>=0 % select the positive component
            Kpos(Kdata(i,1)+1,Kdata(i,2)+1)=Kdata(i,4);
%             Kindspos(Kdata(i,1)+1,Kdata(i,2)+1)=Kdata(i,3);
        end
    end
%      Kneg(logical(isnan(Kneg) .* ~isnan(Kpos)))=0;%-Kpos(boolean(isnan(Kneg) .* ~isnan(Kpos)));
%      Kpos(logical(isnan(Kpos) .* ~isnan(Kneg)))=0;
     Kmat=Kneg; % we like negative component
     Kmat(abs(Kneg+Kpos)>0.005)=nan;
     Kmat(Kmat>0.04)=nan;
     Kmat(Kmat<-0.04)=nan;
%     
    Khist=Kmat(~isnan(Kmat));
%      figure;hist(Khist,40); title ([savedir ' (tn=' num2str(tn) ')'])
     Mtn_meanK=mean(abs(Khist))
%        save ([savedir 'Khist_' num2str(tn) '.mat'],'Khist')
% end

%%
kkk=zeros(mboxsize,nboxsize)*nan;
for i=1:length(Kdata)
    %   kkk(Kdata(i,1)+1,Kdata(i,2)+1)=kkk(Kdata(i,1)+1,Kdata(i,2)+1)+1;
    kkk(Kdata(i,1)+1,Kdata(i,2)+1)=abs(Kdata(i,4));
end
kkk(kkk>0.04)=0.04;

%%
m=2;
kappa=3;
intenergy=1/3*sqrt(2*m*kappa)
Pz=0.020;
hold on
plot([ Pz/intenergy Pz/intenergy], [0 1000],'r','LineWidth',2)
 %%
 figure;
surf(Kmat(1:500,1:500));view([0 90]);
colormap jet;
shading flat
axis equal;
box on;grid off
colorbar;
title ([savedir ' (tn=' num2str(tn) ')'])

%%
%% curvature data for a special grain number
clear
savedir='/home/cenna/Results/2D/Fric800s2_m1_k4/0/';
% start=100;
% step=5000;
% ending=120000;
tn=16000
mboxsize=2000;
nboxsize=mboxsize;
m=1;
kappa=4;
a=sqrt(m/2/kappa);
delx=1;
delt=0.1;
tni=0;
% for tn=[start:step:ending]
%     tni=tni+1;
Kdata=importdata([savedir 'K_' num2str(tn) '.txt']);
grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
gis=find(grains>100); % gis is grain ids that are big
%%
gi=gis(14);
ii=find(Kdata(:,3)==gi);
Kmatgi=zeros(mboxsize,nboxsize)*nan;
for i=1:length(ii)
    Kmatgi(Kdata(ii(i),2)+1,Kdata(ii(i),1)+1)=Kdata(ii(i),4);
end
%  Kmatgi(Kmatgi>0.06)=nan;
%  Kmatgi(Kmatgi<-0.06)=nan;
%
 [in,jn]=find(~isnan(Kmatgi));
        jmin=min(jn);jmax=max(jn);imin=min(in);imax=max(in);
        Kmatgi=Kmatgi(imin:imax,jmin:jmax);
% figure
surf(Kmatgi);view([0 90]);colormap jet;shading flat; axis equal;box on; colorbar;
set(gca,'Clim',[-0.02 0.02]); xlabel('X'); ylabel('Y')
title ([savedir ' (tn=' num2str(tn) ')'])
%%
clear
Msavedir={ '/home/cenna/Results/2D/Fric1000s2_m2_k3_phi_init1000/50/',...
     '/media/Disk2/Results/2D/Fric2000_m2_k3_init2000/30/',...
    '/home/cenna/Results/2D/Fric1000s2_m2_k3_phi_init1000/20_2/',...
           '/media/09004e3f-3e8e-409e-b718-5a0ca5495abd/home/cenna/Results/2D/Fric1200s2_m2_k3_phi/10/',...
    '/home/magnetadmin/Documents/Results/2D/Fric2000_Pz0_m2_k3/',...
    '/home/magnetadmin/Documents/Results/2D/Fric2000_Pz0_m1_k2_run2/',...
    '/home/cenna/Results/2D/Fric2000_m2_k3_init5000/20/',...
    '/home/cenna/Results/2D/Fric800s2_m1_k4/0/',...
    };
Mstart=[500 10000 2000 500 5000 5000 10000 2000];
Mstep=[1000 10000 2000 2000 5000 5000 10000 2000];
Mend=[31500 300000 254000 598500 65000 60000 295000 98000];
MPz=[0.050 0.030 0.020 0.010 0 0 0 0];
figure
for simi=[4]
    Rbar=[];Mtn_meanK=[];
    savedir=Msavedir{simi};
    start=Mstart(simi);
    step=Mstep(simi);
    ending=Mend(simi);
    %     tn=40000
    L=1;
    m=1;
    kappa=2;
    mobility=3/2*L*sqrt(2*kappa/m)
    intenergy=1/3*sqrt(2*m*kappa)
    Pz=MPz(simi);
    delx=1;
    delt=0.1;
    tni=0;
    
    for tn=[start:step:ending]
        tni=tni+1;
        load ([savedir 'Khist_' num2str(tn) '.mat']);
Khist(abs(Khist)>0.05)=[];
%         Khist=HistCurvature(savedir,m,kappa,Pz,tn,'nosave')
         Mtn_meanK(tni)=mean(abs(Khist));
    end

%%
Mtn=[start:step:tn];
timevec=Mtn*delt;
tau=timevec*intenergy*mobility;
% figure
% plot(tau,(1./Mtn_meanK).^2,'o','LineWidth',1.5)
% % plot(timevec,Mtn_meanR.^2,'rs')
% title(savedir)
% xlabel('\tau','FontSize',14);grid on;
% ylabel('(1/<K>)^2','FontSize',14)

%% load grain stats
for tni=1:size(Mtn,2)
    tn=Mtn(tni);
    grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    grains(grains<5)=[];
    Rgrains=sqrt(grains/pi);
    Rbar(tni)=mean(Rgrains)*delx;
    numgrains(tni)=length(grains);
    biggrain(tni)=grains(1)*delx^2;
end
Mtn=[start:step:tn];
timevecG=Mtn*delt;

 figure
% % plotyy(timevec,(1./Mtn_meanK), timevecG, Rbar)
%  plot(1./Mtn_meanK,Rbar(1:length(Mtn_meanK)),':.'); xlabel('1/<K>'); ylabel('<R>')
  plot(1./Rbar(1:length(Mtn_meanK)),Mtn_meanK,':.'); ylabel('<K>'); xlabel('1/<R>')
%% comparison between grain stat and curvature (alpha stuff)
timevec=timevec;
alpha=interp1(timevecG, Rbar, timevec).*Mtn_meanK;
figure
hold on
% plot(Rbar,alpha,'s'); ylabel( '\alpha (<R> \times <K>)','FontSize',14); xlabel(' <R>','FontSize',14)
plot(timevec,alpha,'s'); ylabel( '\alpha (<R> \times <K>)','FontSize',14); xlabel(' Time','FontSize',14)

Ri=linspace(min(Rbar),max(Rbar),10);
alphai=interp1(Rbar,alpha,Ri);
pp=polyfit(Ri,alphai,1)
end
%%
beta=Mtn_meanK*intenergy/Pz

