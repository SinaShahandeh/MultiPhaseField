
%% Hillert Analysis

clear
figure
delx=1;
delt=0.05;
savedir='/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz0_run2/';
L=1;
m=2;
kappa=4;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
start=4000
steps=2000
statstep=60
Mtn=[start:steps:900000];
 Mtn=10000
for tni=1:size(Mtn,2)
    tn=Mtn(tni);
    tnp1=Mtn(tni)+statstep;
    tnm1=Mtn(tni)-statstep;
    grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    
    grainsp1=importdata([savedir 'GrainStat_' num2str(tnp1) '.txt']);
    grainsm1=importdata([savedir 'GrainStat_' num2str(tnm1) '.txt']);
    nonzerogi=find(grains>10);
    grainsvol=grains(nonzerogi)*delx^2;
    grainsvolp1=grainsp1(nonzerogi)*delx^2;
    grainsvolm1=grainsm1(nonzerogi)*delx^2;
    RA=(1/pi)^(1/2)*grainsvol.^(1/2);
    RAp1=(1/pi)^(1/2)*grainsvolp1.^(1/2);
    RAm1=(1/pi)^(1/2)*grainsvolm1.^(1/2);

    GrowthRate=(RAp1-RA)/(1*statstep*delt);
    GrowthRateBar=mean(GrowthRate.*RA)
    Rbar=mean(RA);

    plot(RA/Rbar,GrowthRate.*RA,'.','MarkerSize',8); hold on
    ylabel('(dR/dt) \times R','FontSize',14)
    xlabel('R/<R>','FontSize',14)
    set(gca,'xlim',[0 2.5])
    set(gca,'ylim',[-3 3])
    pp=polyfit(RA./Rbar,GrowthRate.*RA,1);
    Ri=linspace(0, max(RA./Rbar),100);
    Ratei=polyval(pp,Ri);
    plot(Ri,Ratei,'r-')

%         plot(1./RA,GrowthRate,'.','MarkerSize',8); hold on
%         ylabel('Growth Rate (dR_A/dt) ','FontSize',14)
%         xlabel('1/R','FontSize',14)
%         set(gca, 'xlim', [0.001 0.1])
%         set(gca,'ylim',[-0.05 0.05])
%         pp=polyfit(1./RA,GrowthRate,1);
%         Ri=linspace(0, max(1./RA),100);
%         Ratei=polyval(pp,Ri);
%         plot(Ri,Ratei,'r-')

    title(['Time Step= ' num2str(tn)],'FontSize',14);grid on
    set(gca,'LineWidth',2)
    set(gca,'FontSize',14)
    pause(2)
    hold off
end


%% analytical Hillert
figure
L=1;
m=2;
kappa=4;
M=3/2*L*sqrt(2*kappa/m);
sigma=1/3*sqrt(2*m*kappa);
alpha=0.42
Pz=0.030;
R0=11.89;
Rlim=alpha*sigma/Pz;

syms Rcr
f=M/4*alpha*sigma/Rcr*(1-Pz*Rcr/alpha/sigma)^2

Ri=linspace(R0,Rlim-eps,100);
for ni=1:length(Ri)
   ti(ni)=eval(int(1/f,R0,Ri(ni)));
end

hold on
plot(ti,Ri,'g','LineWidth',1.5)
xlabel('Time','FontSize',14);ylabel('<R>', 'FontSize',14)


%% Hillert Analysis with reading all G data  and interpolation of rate

clear
delx=1;
delt=0.1;
savedir='/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz10/'
L=1;
m=2;
kappa=4;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
start=2000
steps=2000
ending=320000
statstep=100
Mtn=[start:steps:ending];
MtnG=[start:statstep:ending];
grains=importdata([savedir 'GrainStat_' num2str(start) '.txt']);
nonzerogi=find(grains>10);
MG=zeros(length(nonzerogi),length(MtnG));
h=waitbar(0,'Reading Grain Stat files...');
for tni=1:length(MtnG)
    tn=MtnG(tni);
    grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    grains=grains(nonzerogi);
    MG(:,tni)=grains;
    waitbar(tni/length(MtnG));
end
close(h)
save([savedir 'AllGrainStat.mat'])
%% 
clear
savedir='/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz10/';
load([savedir 'AllGrainStat.mat']);
%% Kinetics and topology of 'a' grain 

initgrains=importdata([savedir 'GrainStat_' num2str(start) '.txt']);
initnonzerogi=find(initgrains>10);

rangestep=200
gii=241 % index in the MG matrix
subplot(3,1,1);
grainsvol=MG(gii,:)*delx^2;
%     R=sqrt(grainsvol/pi);   %%%%%%%%%%%% 2D %%%%%%%%
     R=(3*grainsvol/4/pi).^(1/3);  %%%%%%%%%%%% 3D %%%%%%%%
plot(MtnG*delt,R);grid on
ylabel('R','FontSize',14)
Mtn=[start+2000:2000:ending-2000];
GrowthRate=zeros(1,length(Mtn));
RMtn=[];GrowthRate=zeros(size(Mtn));FaceMtn=zeros(size(Mtn));
for tni=1:size(Mtn,2)
    tn=Mtn(tni);
    MGtnind=find(MtnG==tn);
    t=[tn-rangestep:statstep:tn+rangestep]*delt;
    Gi=MG(gii,MGtnind-rangestep/statstep:MGtnind+rangestep/statstep);
%     Gi=sqrt(Gi/pi);        %%%%%%%%%%%% 2D %%%%%%%%
    Gi=(3*Gi/4/pi).^(1/3); %%%%%%%%%%%% 3D %%%%%%%%
    pp=polyfit(t,Gi,2);
    GrowthRate(tni)=polyval(polyder(pp),tn*delt)/2;
    RMtn(tni)=sqrt(MG(gii,MGtnind)*delx^2/pi);
    
    load([savedir 'facenums/facenum_' num2str(tn) '.mat']);
    nonzerofacegi=find(MG(:,MGtnind)>10);
    indface=find(nonzerofacegi==gii);
    if isempty(indface)==true
        FaceMtn(tni)=0;
    else
        FaceMtn(tni)=facenum(indface);
    end
end
subplot(3,1,2)
plot(Mtn*delt,GrowthRate.*RMtn/mobility/intenergy,'x','MarkerSize',6,'LineWidth',1); grid on
ylabel('(dR/dt) \times R/M / \sigma','FontSize',14)
subplot(3,1,3)
plot(Mtn*delt,FaceMtn,'^','MarkerSize',6,'LineWidth',1); grid on
ylabel('N','FontSize',14)
xlabel('Time','FontSize',14)

%% kinetics map R vs. dR/dt
initgrains=importdata([savedir 'GrainStat_' num2str(start) '.txt']);
initnonzerogi=find(initgrains>10);

rangestep=200
gii=110 % index in the MG matrix
Mtn=[start+200:100:ending-200];
GrowthRate=zeros(1,length(Mtn));
RMtn=[];GrowthRateMtn=zeros(size(Mtn));
for tni=1:size(Mtn,2)
    tn=Mtn(tni);
    MGtnind=find(MtnG==tn);
    t=[tn-rangestep:statstep:tn+rangestep]*delt;
    Gi=MG(gii,MGtnind-rangestep/statstep:MGtnind+rangestep/statstep);
    Gi=sqrt(Gi/pi);        %%%%%%%%%%%% 2D %%%%%%%%
    %         Gi=(3*Gi/4/pi).^(1/3); %%%%%%%%%%%% 3D %%%%%%%%
    pp=polyfit(t,Gi,2);
    GrowthRateMtn(tni)=polyval(polyder(pp),tn*delt);
    RMtn(tni)=sqrt(MG(gii,MGtnind)*delx^2/pi);
    nonzerogi=find(MG(:,MGtnind)>10);
    grainsvol=MG(nonzerogi,MGtnind);
        RNZ=sqrt(grainsvol/pi);   %%%%%%%%%%%% 2D %%%%%%%%
%   RNZ=(3*grainsvol/4/pi).^(1/3);  %%%%%%%%%%%% 3D %%%%%%%%
    Rbar(tni)=mean(RNZ);
end
figure
% plot(RMtn./Rbar,GrowthRateMtn.*RMtn/mobility/intenergy,'.','MarkerSize',6,'LineWidth',1); grid on
comet(RMtn,GrowthRateMtn.*RMtn/mobility/intenergy,0.001)
% ylabel('(dR/dt) \times R/M / \sigma','FontSize',14)
% xlabel('R/<R>','FontSize',14)
%% Hillert Scatter
clear
savedir='/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.010_m1_k2_init4000_run1/';
load([savedir 'AllGrainStat.mat']);
figure
Mtn=[start+200:steps:ending];
Mtn=4100
for tni=1:size(Mtn,2)
    tn=Mtn(tni);
    MGtnind=find(MtnG==tn);
    grains=MG(:,MGtnind);
%     nonzerogi=find(grains>10);
    nonzerogi=1:length(grains); % select all grains
    rangestep=100
    GrowthRate=zeros(length(nonzerogi),1);
    for gi=1:length(nonzerogi)
        gii=nonzerogi(gi);
        t=[tn-rangestep:statstep:tn+rangestep]*delt;
        Gi=MG(gii,MGtnind-rangestep/statstep:MGtnind+rangestep/statstep);
%         Gi=sqrt(Gi/pi);        %%%%%%%%%%%% 2D %%%%%%%%
         Gi=(3*Gi/4/pi).^(1/3); %%%%%%%%%%%% 3D %%%%%%%%
       pp=polyfit(t,Gi,2);
       GrowthRate(gi)=polyval(polyder(pp),tn*delt);
%        figure;plot(t,Gi,'-*')
%        GrowthRate(gi)*R(gi)
    end
    nonzerogi=find(MG(:,MGtnind)>10);
    GrowthRateNZ=GrowthRate(nonzerogi);
    GrowthRateBar=mean(GrowthRateNZ);
    grainsvol=MG(nonzerogi,MGtnind);
%    R=sqrt(MG(:,MGtnind)/pi);   %%%%%%%%%%%% 2D %%%%%%%%
    R=(3*MG(:,MGtnind)/4/pi).^(1/3);  %%%%%%%%%%%% 3D %%%%%%%%
%    RNZ=sqrt(grainsvol/pi);   %%%%%%%%%%%% 2D %%%%%%%%
    RNZ=(3*grainsvol/4/pi).^(1/3);  %%%%%%%%%%%% 3D %%%%%%%%
    Rbar=mean(RNZ);
    plot(RNZ/Rbar,GrowthRateNZ.*RNZ/mobility/intenergy,'.','MarkerSize',8); hold on
    ylabel('R (dR/dt)/M / \sigma','FontSize',14)
    xlabel('R/<R>','FontSize',14)
    set(gca,'xlim',[0 2.5])
%    set(gca,'ylim',[-0.5 0.5]) %%%%%%%%%%%% 2D %%%%%%%%
    set(gca,'ylim',[-0.9 0.9]) %%%%%%%%%%%% 3D %%%%%%%%
    
%     pp=polyfit(RNZ./Rbar,GrowthRateNZ.*RNZ/mobility/intenergy,1);
%     Ri=linspace(0, max(R./Rbar),100);
%     Ratei=polyval(pp,Ri);
%     plot(Ri,Ratei,'r-')

    
    title(['Time Step= ' num2str(tn)],'FontSize',14);grid on
    set(gca,'LineWidth',2)
    set(gca,'FontSize',14)
    pause(0.5)
    hold off
end

%% Calculate Rcr with from grain size distribution
Pz=0.010;
alpha=1.2;
z=Pz/intenergy;

udata=RNZ/Rbar;
u=linspace(0.01,2.5,500); 
% from wbl fit
wbl_p=wblfit(udata);
P=pdf('wbl',u,wbl_p(1),wbl_p(2));
Ri=u.*Rbar; % because I got P from u that was obtained by deviding RNZ to Rbar
Rcr=fsolve(@(Rcr) integral(Ri,P,Rcr,z,alpha),10)
Rcr/Rbar
%% analytical hillert with pinning
K=alpha*mobility*intenergy;
Ri=linspace(min(R),5*Rbar,300);
    for i=1:length(Ri)

        dRdt(i)=0; % else when pinning is stronger
         if 1/Ri(i)>=1/Rcr+z/alpha
             dRdt(i)=K*(1/Rcr-1/Ri(i)+z/alpha);
         end
         if 1/Ri(i)<1/Rcr-z/alpha
             dRdt(i)=K*(1/Rcr-1/Ri(i)-z/alpha);
         end
    end

hold on
plot(Ri/Rbar,dRdt.*Ri/intenergy/mobility,'k','LineWidth',1.5)


[AX,H1,H2]=plotyy(Ri/Rbar,dRdt.*Ri/intenergy/mobility,u,P)
set(gca,'YTick',[-1:0.5:1])
set(H1,'LineWidth',1.5,'Color','r')
set(H1,'LineWidth',1,'Color','r')
title(['Time Step= ' num2str(tn) ', R_{cr}=' num2str(Rcr/Rbar)],'FontSize',14);grid on
%% colour coded Hillert
clear
savedir='/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.010_m1_k2_init4000_run1/';
load([savedir 'AllGrainStat.mat']);
L=1;
m=1;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
figure
Mtn=[start+200:steps:ending];
Mtn=6000
for tni=1:size(Mtn,2)
    tn=Mtn(tni);
    MGtnind=find(MtnG==tn);
    grains=MG(:,MGtnind);
%     nonzerogi=find(grains>10);
    nonzerogi=1:length(grains); % select all grains
    rangestep=100
    GrowthRate=zeros(length(nonzerogi),1);
    for gi=1:length(nonzerogi)
        gii=nonzerogi(gi);
        t=[tn-rangestep:statstep:tn+rangestep]*delt;
        Gi=MG(gii,MGtnind-rangestep/statstep:MGtnind+rangestep/statstep);
%         Gi=sqrt(Gi/pi);        %%%%%%%%%%%% 2D %%%%%%%%
                 Gi=(3*Gi/4/pi).^(1/3); %%%%%%%%%%%% 3D %%%%%%%%
        pp=polyfit(t,Gi,2);
        GrowthRate(gi)=polyval(polyder(pp),tn*delt);
        %        figure;plot(t,Gi,'-*')
        %        GrowthRate(gi)*R(gi)
    end
    nonzerogi=find(MG(:,MGtnind)>10);
    GrowthRateNZ=GrowthRate(nonzerogi);
    GrowthRateBar=mean(GrowthRateNZ);
    grainsvol=MG(nonzerogi,MGtnind);
%     R=sqrt(MG(:,MGtnind)/pi);   %%%%%%%%%%%% 2D %%%%%%%%
    R=(3*MG(:,MGtnind)/4/pi).^(1/3);  %%%%%%%%%%%% 3D %%%%%%%%
%     RNZ=sqrt(grainsvol/pi);   %%%%%%%%%%%% 2D %%%%%%%%
    RNZ=(3*grainsvol/4/pi).^(1/3);  %%%%%%%%%%%% 3D %%%%%%%%
    Rbar=mean(RNZ);

    pp=polyfit(RNZ./Rbar,GrowthRateNZ.*RNZ/mobility/intenergy,1);
    Ri=linspace(0, max(R./Rbar),100);
    Ratei=polyval(pp,Ri); % plot an interpolaton line for Hillert fit.
%     plot(Ri,Ratei,'k-','LineWidth',1.5); hold on
    % facenum loads face num that matches the R above for grains larger than 10.
    load([savedir 'facenums/facenum_' num2str(tn) '.mat']);

%     rgb=[0 0 0 0.1 0.5 0.6 0.7 0.3 0.1 0.9 0.9 0.2; 0 0 0 0.1 0.6 0.4 0.3 0.2 0.5 0.8 0 0.99; 0 0 0 0.99 0.5 0.8 0.6 0.1 0.5 0.8 0 0.2];
         rgb=[rand(1,30); rand(1,30); rand(1,30)];
    % plotting Hillert Scatter
%     for n=3:30
%         indn=find(facenum==n);
%         RClass=RNZ(indn);
%         dRdtClass=GrowthRateNZ(indn);
%         plot(RClass/Rbar,dRdtClass.*RClass/mobility/intenergy,'.','MarkerSize',8,'color',[rgb(1,n) rgb(2,n) rgb(3,n)]); hold on
%     end
%     ylabel('(dR/dt) \times R/M / \sigma','FontSize',14)
%     xlabel('R/<R>','FontSize',14)
%     set(gca,'xlim',[0 2.5])
%     set(gca,'ylim',[-1 1])
%     plot([0 2.5],[-1 1*(8/9*2.5-1)],'r','LineWidth',1.5)
    % plotting von Neumann scatter
    i=1;
%     Ni=3:12; %%%%%%%%%%%% 2D %%%%%%%%
    Ni=4:30; %%%%%%%%%%%% 3D %%%%%%%%
    for n=Ni
        indn=find(facenum==n);
        RClass=RNZ(indn);
        meanRClass(i)=mean(RClass);
        dRdtClass=GrowthRateNZ(indn);
        plot(zeros(size(RClass))+n,dRdtClass.*RClass/mobility/intenergy,'.',...
            'MarkerSize',8,'color',[rgb(1,n) rgb(2,n) rgb(3,n)]); hold on
        vonY(i)=mean(dRdtClass.*RClass/mobility/intenergy);
        i=i+1;
    end
    plot(Ni,vonY,'s','LineWidth',2)
    ylabel('R (dR/dt)/M / \sigma','FontSize',14)
    xlabel('N','FontSize',14)
    set(gca,'xlim',[3 12])
    set(gca,'ylim',[-0.5 0.5])
%     plot([0 12],[-1 1],'k','LineWidth',1.5)
    box on
    set(gca,'LineWidth',2)
    set(gca,'FontSize',14)
    pause(0.5)
    grid on
    title(['Time step='  num2str(tn)])
end
%% R - N
figure
plot(facenum,RNZ/Rbar,'.')
hold on

plot([3:12],meanRClass/Rbar,'ro')
pp=polyfit([3:10],meanRClass(1:8)/Rbar,2)
%% analytical von-neumann with pinning
Pz=0.01;
z=(7.15/2/pi)*6*Pz/intenergy;
K=mobility*intenergy/6;
Rcr=1*Rbar
Ni=3:12;
% Ri=Rcr*(Ni-6)/6/0.5+Rcr;
Ri=Rbar*(0.2305*Ni-0.4541);
% Ri=Rbar*polyval(pp,Ni);
    for ni=1:length(Ni)
        RdRdt(ni)=0; % else when pinning is stronger
         if Ni(ni)<6-z*Ri(ni)
             RdRdt(ni)=K*(Ni(ni)-6+z*Ri(ni));
         end
         if Ni(ni)>6+z*Ri(ni)
             RdRdt(ni)=K*(Ni(ni)-6-z*Ri(ni));
         end
    end

hold on
plot(Ni,RdRdt/intenergy/mobility,'k','LineWidth',1.5)

