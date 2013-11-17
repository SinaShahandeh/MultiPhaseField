clear
savedir='/home/magnetadmin/Documents/Results/3D/Fric300_Pz0_m1_k2_init0/';
start=500;
step=500;
ending=350000;

mboxsize=300;
nboxsize=mboxsize;
kboxsize=mboxsize;

tn=15000
tni=0;
% for tn=[start:step:ending]
%     tni=tni+1;
    Kdata=importdata([savedir 'K_' num2str(tn) '.txt']);
    Kneg=zeros(mboxsize,nboxsize,nboxsize)*nan;
%     Kinds=zeros(mboxsize,nboxsize,nboxsize)*nan;
    for i=1:length(Kdata)
        if Kdata(i,5)<=0 % select the negative component
            Kneg(Kdata(i,1)+1,Kdata(i,2)+1,Kdata(i,3)+1)=Kdata(i,5);
%             Kinds(Kdata(i,1)+1,Kdata(i,2)+1,Kdata(i,3)+1)=Kdata(i,4);
        end
    end
    Kpos=zeros(mboxsize,nboxsize,nboxsize)*nan;
    for i=1:length(Kdata)
        if Kdata(i,5)>=0 % select the positive component
            Kpos(Kdata(i,1)+1,Kdata(i,2)+1,Kdata(i,3)+1)=Kdata(i,5);
        end
    end
%      Kneg(logical(isnan(Kneg) .* ~isnan(Kpos)))=0;%-Kpos(boolean(isnan(Kneg) .* ~isnan(Kpos)));
%      Kpos(logical(isnan(Kpos) .* ~isnan(Kneg)))=0;
     Kmat=Kneg; % we like negative component
     Kmat(abs(Kneg+Kpos)>0.002)=nan;
     Kmat(Kmat>0.04)=nan;
     Kmat(Kmat<-0.04)=nan;
%     
    Khist=Kmat(~isnan(Kmat));
    Mtn_meanK=mean(abs(Khist))
%     figure;hist(abs(Khist/Mtn_meanK),40); title ([savedir ' (tn=' num2str(tn) ')'])
%     x=linspace(0,0.035/Mtn_meanK,40); % bins
%     [N,X]=hist (abs(Khist)/Mtn_meanK,x); %
%     plot(X,N/length(Khist),'LineWidth',2);
%     hold on
    
%    save ([savedir 'Khist_' num2str(tn) '.mat'],'Khist')
% end
% Kpos(Kpos>0.01)=nan;
% sliceomatic(Kneg)
%% %% cross section on curvature map
figure;
slice(Kpos,[1 ],[1],[kboxsize ]);
% axis([1 (mboxsize) 1 nboxsize 1 kboxsize]);
% surf(abs(Kpos(:,:,5)));view([0 90]);
colormap jet; shading flat;
axis equal; box off; grid off
colorbar; set(gca,'Xlim',[0 length(Kmat)]); set(gca,'Ylim',[0 length(Kmat)])
set(gca,'LineWidth',1.5),set(gca,'Clim',[0 0.03])
% title ([savedir ' (tn=' num2str(tn) ')'])
% mkdir([savedir 'curvarure_map/'])
% print([savedir 'curvarure_map/' num2str(tn) '.png'],'-dpng','-r300',gcf)
sdsd
    %% Curvature Histograms
savedir='/home/magnetadmin/Documents/Results/3D/Fric512_Pz0_m1_k2_init0/';
tn=10000
m=1
kappa=2
Pz=0
maxK=0.025;
dosave='nosave';
PlotType='plot';
figure
Khist=HistCurvature(savedir,m,kappa,Pz,tn,maxK,dosave,PlotType);
 %% Raw curvature histograms
clear
savedir='/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.012_m1_k2_init4000/';
tn=900
Pz=0.015
L=1;
m=1;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m)
intenergy=1/3*sqrt(2*m*kappa)
tni=0;
start=60050;
step=500;
ending=100000000;
figure
 for tn=[start:step:ending]
    tni=tni+1;
    Kdata=importdata([savedir 'K_' num2str(tn) '.txt']);
    Khist=Kdata(:,5);
    Khist(Khist>0)=[]; %----> Selection Criteria
    Khist(Khist<-0.1)=[];
    Khist=abs(Khist);
    meanK(tni)=mean(Khist);
    save( [savedir 'Khist_neg_' num2str(tn) '.mat'], 'Khist')
end
%%
Mi=Kdata(:,1);
Mj=Kdata(:,2);
Mk=Kdata(:,3);
ni=0;
for ki=1:length(Mi)
    reppoints=find([(Mi==Mi(ki)) & (Mj==Mj(ki)) & (Mk==Mk(ki))]==1);
    if length(reppoints)==4 % that is a grain boundary face
  ni=ni+1;
        
    end
end

%% curvature data for a special grain number
clear
savedir='/home/magnetadmin/Documents/Results/3D/Fric100_Pz0.010_m1_k2_init1000/';
start=5000;
step=5000;
 ending=130000;
tn=5000
mboxsize=100;
nboxsize=mboxsize;
kboxsize=mboxsize;
m=1;
kappa=2;
a=sqrt(m/2/kappa);
delx=1;
delt=0.1;
tni=0;
% for tn=[start:step:ending]
%     tni=tni+1;
Kdata=importdata([savedir 'K_' num2str(tn) '.txt']);
grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
gis=find(grains>10); % gis is grain ids that are big
gi=gis(1);
ii=find(Kdata(:,4)==gi);
Kmatgi=zeros(mboxsize,nboxsize,nboxsize)*nan;
for i=1:length(ii)
    Kmatgi(Kdata(ii(i),1)+1,Kdata(ii(i),2)+1,Kdata(ii(i),3)+1)=Kdata(ii(i),5);
end
Kmatgi(Kmatgi>0.06)=nan;
sliceomatic(Kmatgi)
%
Khisti=Kmatgi(~isnan(Kmatgi));
figure
hist(Khisti,50);
title(savedir)
title (['Total Sample points = ' num2str(length(Khisti)) ])
AverageCurvature=mean(Khisti)
%% Statistics :sampling based on the area of the segments
%% Curvature of a segment
Lseg=bwlabeln(~isnan(Kmat));
i=9
Segi=find(Lseg==i);
Segimat=Kmat.*[Lseg==i];
sliceomatic(Segimat)
[f,ki]=hist((Kmat(Segi)),30); % finding position of maximum in histogram distribution and repeating the data for the whole segment
figure; bar(ki,f)
[a, fmaxi]=max(f)
(ki(fmaxi))
%% Total curvature of the structure
meanK=[];
Khist=[];
Lseg=bwlabeln(~isnan(Kmat));
for i=1:max(max(max(Lseg)))
    Segi=find(Lseg==i);
    if length(Segi)>100 % removing very small segments
        Ksegmat=Kmat(Segi);
         [f,ki]=hist(Ksegmat,500); % finding position of maximum in histogram distribution and repeating the data for the whole segment
%           hist(log(Kmat(Segi)),100); title([ num2str(i) ' -- <K>= ' num2str(mean(Kmat(Segi))) ' -- A= ' num2str(length(Segi)) ] ); pause (2);
         [a, fmaxi]=max(f);
%         Ksegi=(ki(fmaxi))+zeros(length(Segi),1);
         rangek=Ksegmat(Ksegmat<ki(fmaxi)+0.01 & Ksegmat>ki(fmaxi)-0.01);
         Ksegi=mean(rangek)+zeros(length(Segi)+fix(2*sqrt(pi)*sqrt(length(Segi)/2)*0^2),1);
        Khist=[Khist Ksegi'];
    end
end
figure
hist((Khist),30)
title ([  'Ttotal Sample points = ' num2str( Khist) ])
meanK=[meanK mean(Khist)]

%% ------------------------------------------------------------------------
%% Sampling of curvature data from Kmat over time steps (kmat_seg_neg GENERATOR)
clear
savedir='/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.015_m1_k2_init2000/';
start=5000;
step=5000;
ending=125000;
mboxsize=300;
nboxsize=300;
kboxsize=300;
m=1;
kappa=2;
a=sqrt(m/2/kappa);
delx=1;
delt=0.1;
tni=0;
for tn=[start:step:ending]
    tni=tni+1;
    Kdata=importdata([savedir 'K_' num2str(tn) '.txt']);
    Kmat=zeros(mboxsize,nboxsize,nboxsize)*nan;
    for i=1:length(Kdata)
        if Kdata(i,5)<0 % -----> select the positive or negative component
            Kmat(Kdata(i,1)+1,Kdata(i,2)+1,Kdata(i,3)+1)=Kdata(i,5);
        end
    end
    Kmat=abs(Kmat);
    Kmat(Kmat>0.1)=nan; % to detach segments at the triple junctions with high curvature
    Khist=[];
    Lseg=bwlabeln(~isnan(Kmat));
    segnum(tni)=0;
    for i=1:max(max(max(Lseg)))
        Segi=find(Lseg==i);
        if length(Segi)>1 % removing very small segments
            Ksegmat=Kmat(Segi);
            segnum(tni)=segnum(tni)+1;
            [f,ki]=hist((Kmat(Segi)),500); % finding position of maximum in histogram distribution and repeating the data for the whole segment
            [a, fmaxi]=max(f);
            rangek=Ksegmat(Ksegmat<ki(fmaxi)+0.01 & Ksegmat>ki(fmaxi)-0.01);
            Ksegi=mean(rangek)+zeros(length(Segi),1);
            Khist=[Khist Ksegi'];
        end
    end
    meanK(tni)=mean(Khist);
    figure
    hist((Khist),100)
    save( [savedir 'Khist_seg_neg' num2str(tn) '.mat'], 'Khist')
end

%% Curvature Histograms
savedir='/home/magnetadmin/Documents/Results/3D/Fric256_Pz0.010_m1_k2_init4000_trans/';
start=500;
step=500;
ending=130000;
Mtn=[start:step:ending];
Mtn=4000
figure
for tn=Mtn
m=1
kappa=2
intenergy=1/3*sqrt(2*m*kappa);
Pz=0.010
maxK=0.030;
maxY=0.15;
dosave='save'; % 'nosave' or 'save'
PlotType='bar'; % bar or plot

hold off
Khist=HistCurvature(savedir,m,kappa,Pz,tn,maxK,maxY,PlotType);
pause(0.1)
if strcmp(dosave,'save')==true;
    mkdir([savedir 'curvarure_dist_Khist/'])
    print([savedir 'curvarure_dist_Khist/' num2str(tn) '.pdf'],'-dpdf','-r200',gcf)
end
end
%% reading Khist mat files and plot distributions

%%
figure
Ktime=[start:step:ending]*delt;
plot(Ktime, (meanK),'ro')
xlabel('Time')
ylabel('(Curvature Based Equivalent Radius)^2,  (1/K)^2')
hold on
plot([Ktime(1) Ktime(end)], (1./[Pz/intenergy Pz/intenergy]).^2,'LineWidth',1.5 , 'Color','r')
%% load grain stats

Mtn=[start:step:ending];
try
    for tni=1:size(Mtn,2)
        tn=Mtn(tni);
        GrainStat=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
        grains=GrainStat;
        grains(grains<3^3)=[];
        radius=((3/4/pi)^(1/3))*grains.^(1/3);
        Rbar(tni)=mean(radius);
        numgrais(tni)=length(grains);
    end
catch
end
Mtn=[start:step:tn];
timevec=Mtn*delt;

figure
plotyy(Ktime,(2./meanK), timevec, Rbar)

%% comparison between grain stat and curvature (alpha stuff)

alpha=interp1(timevec, Rbar, Ktime,'spline').*meanK
figure
% plot(Ktime,alpha,'s')
plot(Rbar,alpha,'d')
ylabel( '\alpha (Mean Curvature \times Equivalent Radius)')
xlabel(' Time')

%% --------------------------- all steps for alpha and beta -----------

clear
Msavedir={ ...
'/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.010_m1_k2_init4000/', ...
'/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.012_m1_k2_init4000/', ...
'/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.013_m1_k2_init4000/', ...
'/home/magnetadmin/Documents/Results/3D/Fric256_Pz0.013_m1_k2_init4000_trans/', ...
'/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.015_m1_k2_init2000/', ...
'/home/magnetadmin/Documents/Results/3D/Fric300_Pz0_m1_k2_init0_run2/', ...
'/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.012_m1_k2_init4000_trans2/', ...
'/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.012_m1_k2_init4000_trans/', ...
'/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.012_m1_k2_cont30000/' , ...
    };
figure
Mstart=[5000 5000 5000 500 5000 5000 500 4000 30000];
Mstep=[5000 5000 5000 500 5000 5000 500 500 2000];
for sim=[ 2]
  savedir=Msavedir{sim};
start=Mstart(sim);
step=Mstep(sim);
ending=1000000;
delt=0.1;
tni=0;
L=1;
m=1;
kappa=2
mobility=3/2*L*sqrt(2*kappa/m)
intenergy=1/3*sqrt(2*m*kappa)
% Pz=0.010
clear Rbar meanK
try
for tn=[start:step:ending]
    tni=tni+1;
    load([savedir 'Khist_neg_' num2str(tn) '.mat'])
    meanK(tni)=mean(Khist);
end
end
ending=tn-step;
Ktime=[start:step:ending]*delt;
% load grain stats
% step=100;
Mtn=[start:step:ending];
try
    for tni=1:size(Mtn,2)
        tn=Mtn(tni);
        GrainStat=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
        grains=GrainStat;
        grains(grains<3^3)=[];
        radius=((3/4/pi)^(1/3))*grains.^(1/3);
        Rbar(tni)=mean(radius);
    end
catch
end
Mtn=[start:step:tn];
timevec=Mtn*delt;
% comparison between grain stat and curvature (alpha stuff)
alpha=interp1(timevec, Rbar, Ktime,'spline').*meanK;

hold on
plot(Rbar,alpha,'d');xlabel(' <R>'); ylabel( '\alpha (<K> \times <R>)')
%  plot(Rbar,1./meanK,'s'); xlabel(' Rbar');ylabel('<K>^{-1}')
% plot(1./meanK,Rbar,'s'); ylabel(' Rbar');xlabel('<K>^{-1}')
% plot(Ktime,meanK,'bo'); xlabel('Time'); ylabel('<K>');grid on
% plot(timevec,Rbar,'d'); xlabel('Time'); ylabel('<R>');grid on
end
% ----------------------------------------------------
%% calculating beta  :  \sigma K = \beta Pz

beta=intenergy*meanK/Pz
figure
plot(Ktime, beta,'*')
beta=beta(end)



%%


%% ------------------------------- Close looks---------------
%% analytical sphere and its curvature
clear
m=1;
kappa=4;
a=sqrt(m/2/kappa);
r0=20;
syms r
e=1/2*(1-tanh(a*(r-r0)));
% ezplot(e, [0 50])
syms x y
etasph=subs(e,'r','sqrt(x^2+y^2+z^2)')
eta_x=diff(etasph,'x');
eta_y=diff(etasph,'y');
eta_z=diff(etasph,'z');
eta_xx=diff(eta_x,'x');
eta_yy=diff(eta_y,'y');
eta_zz=diff(eta_z,'z');
gradeta=-0.5*a*(1-(2*etasph-1)^2);
grad2eta=-a^2*(2*etasph-1)*(1-(2*etasph-1)^2);
del2eta=eta_xx+eta_yy+eta_zz;

H2=simple((del2eta-grad2eta)/gradeta)
% ezmesh(H,[-50 50],[-50 50])
%%
clear
m=1;
kappa=4;
a=sqrt(m/2/kappa);
r0=20;
x=linspace(-35,35,70);
y=x;
z=x;
[X,Y,Z]=meshgrid(x,y,z);
num_eta=zeros(length(x),length(y),length(z));
num_eta=1/2*(1-tanh(a*(sqrt(X.^2+Y.^2+Z.^2)-20)));

% sliceomatic(num_eta,x,y,z)
etani=num_eta;
%% calculating curvature
L=1;
m=1;
kappa=4;
mobility=3/2*L*sqrt(2*kappa/m)
intenergy=1/3*sqrt(2*m*kappa)
mboxsize=length(x);nboxsize=mboxsize;kboxsize=nboxsize;
delx=(max(x)-min(x))/length(x);
for i=1:mboxsize
    for j=1:nboxsize
        for k=1:kboxsize
            sumeta=etani(indg(i+1,mboxsize),j,k)+etani(indg(i-1,mboxsize),j,k)+etani(i,indg(j+1,nboxsize),k)+etani(i,indg(j-1,nboxsize),k)...
                +etani(i,j,indg(k+1,kboxsize))+etani(i,j,indg(k-1,kboxsize));
             sumeta2=etani(indg(i+1,mboxsize),indg(j+1,nboxsize),k)+etani(indg(i-1,mboxsize),indg(j+1,nboxsize),k)+etani(indg(i+1,mboxsize),indg(j-1,nboxsize),k)+etani(indg(i-1,mboxsize),indg(j-1,nboxsize),k)...
                 +etani(i,indg(j+1,nboxsize),indg(k+1,kboxsize))+etani(i,indg(j-1,nboxsize),indg(k+1,kboxsize))+etani(i,indg(j+1,nboxsize),indg(k-1,kboxsize))+etani(i,indg(j-1,nboxsize),indg(k-1,kboxsize))...
                 +etani(indg(i+1,mboxsize),j,indg(k+1,kboxsize))+etani(indg(i-1,mboxsize),j,indg(k+1,kboxsize))+etani(indg(i+1,mboxsize),j,indg(k-1,kboxsize))+etani(indg(i-1,mboxsize),j,indg(k-1,kboxsize));
             del2eta(i,j,k)=1/6/(delx*delx)*(sumeta2+2*sumeta-24*etani(i,j,k));
%             del2eta(i,j,k)=1/(delx*delx)*(sumeta-6*etani(i,j,k));
        end
    end
end
%%

gradeta=-0.5*sqrt(m/2/kappa)*(1-(2*etani-1).^2);
grad2eta=-sqrt(m/2/kappa)^2*(2*etani-1).*(1-(2*etani-1).^2);
H=(del2eta-grad2eta)./gradeta;
% H(abs(H)>0.5)=nan;
H(etani>0.8)=nan;
H(etani<0.2)=nan;
sliceomatic(H,x,y,z)
figure
midind=length(x)/2;
plot(x,del2eta(:,midind,midind))
hold on
plot(x,gradeta(:,midind,midind),'r')
plot(x,grad2eta(:,midind,midind),'k')
plot(x,H(:,midind,midind),'b')
plot(x,etani(:,midind,midind),'g:')
grid; xlabel('X'); ylabel('Profile')
figure
hsamples=H(~isnan(H));
hist(hsamples,50)

%%
figure
subplot(2,1,1)
plot(eta(:,nboxsize/2),'.-')
hold on
xi=[0:mboxsize-1]*delx; x1=35;
analyticeta=0.5*(1+tanh(sqrt(m/2/kappa)*(x1-xi)));
plot(analyticeta,'r.-')
title(num2str(tn));ylabel('\eta')
hold off
subplot(2,1,2)
plot(H(:,nboxsize/2),'.-'); hold on;
hold off; ylabel('Local Curvature'); xlabel('Position');grid on; axis([0 mboxsize -0.05 -0.03])

%%

%% statistic sampling of curvature map

[i,j]=find(~isnan(H));
for n=1:length(i)
    Hsamples(n)=H(i(n),j(n));
end
mean(Hsamples)
figure
hist(Hsamples)
xlabel('Curvature Data')
ylabel('Counut')
