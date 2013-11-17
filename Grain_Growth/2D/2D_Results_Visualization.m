%% ------------------------- 2D Grain Growth ----------------------------

%% 2D cross sections

clear
ni=0
figure
tn=30000;
%  for tn=[1000:1000:600000]
ni=ni+1;
savedir='/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz10_run2/';
phi=importdata([savedir 'Fullres_' num2str(tn) '.txt']);
%  figure
% imshow(phi);
surf(phi);
view ([0 90])
shading flat
axis off
box off
caxis([0 1]);
colormap gray
shading flat
axis off
box off
title(['Time Step= ' num2str(tn)])
set(gca,'DataAspectratio',[1 1 1])
pause(0.01);
%    mkdir([savedir 'phi/'])
%    print([savedir 'phi/' num2str(ni) '.png'],'-dpng','-r200',gcf)
%    imwrite(phi,[savedir '/' num2str(ni) '.png'],'png')
%   end

%% Reading Inds
clear
figure
tn=2000;
% for tn=[10:10:20000000]
savedir='/media/Disk2/Results/2D/Fric2000_m2_k6_init2000/15/';
phi=importdata([savedir '/Inds_' num2str(tn) '.txt']);
surf(phi);
view ([0 90])
shading flat
axis off
box off
title(['Time= ' num2str(tn)])
set(gca,'DataAspectratio',[1 1 1])
pause(0.05);
colorbar
% end

%% Recrytallization structure
clear
figure
start=500;
step=500;
ending=500
ni=0;
for tn=[start:step:ending]
    ni=ni+1;
delt=0.15;
savedir='/home/cenna/Results/2D/flat_grain/0/';
phi=importdata([savedir '/Fullres_' num2str(tn) '.txt']);
inds=importdata([savedir '/Inds_' num2str(tn) '.txt']);
im(:,:,1)=phi;
im(:,:,2)=phi.*[inds~=1];
im(:,:,3)=phi.*[inds~=1];
imshow(im);
caxis([0 1]);
colormap gray
shading flat
axis off
box off
title(['Time= ' num2str(tn)])
set(gca,'DataAspectratio',[1 1 1])
pause(01);
print([savedir '/' num2str(ni) '.png'],'-dpng','-r200',gcf)
end

%% recrystallization fraction
clear
mboxsize=2000;
nboxsize=mboxsize;
delx=1;

savedir='/home/cenna/Results/2000/';
data=importdata([savedir '/MatrixA' '.txt']);
time=data(:,1);
matrixa=data(:,2);
frac=1-(matrixa/(mboxsize*nboxsize*delx^2));
figure
plot(time,frac)

incubationtime=50
X=log(time(frac<0.95)+incubationtime);
frac=frac(frac<0.95);
Y=log(log(1./(1-frac)));
figure
plot(X,Y)


%% Statistics from GrainStats
clear
Msavedir={'/media/09004e3f-3e8e-409e-b718-5a0ca5495abd/home/cenna/Results/2D/Fric1200s2_m2_k3_phi/10/', ...
    '/media/Disk2/Results/2D/Fric2000_Pz0_m2_k3/', ...
           '/home/cenna/Results/2D/Fric1000s2_m2_k3_phi_init1000/20_2/',...
           '/media/Disk2/Results/2D/Fric2000_Pz0_m1_k2_run2/' , ...
    '/media/Disk2/Results/2D/solute_m2_k4/Fric2000_a0_b1000_m2_k4/',...
    '/home/cenna/Results/2D/Fric800s2_m1_k4/0/',...
    '/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz0/',...
    '/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz10/',...
    '/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz20/',...
    '/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz30/',...
    '/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz10_run2/',...
    };
Mstart=[ 1000 2000 2000 2000 2000 4000 2000 2000 2000 2000 2000 2000 2000];
Msteps=[6000 1000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000];
Mend=[600000 ones(1,11)*600000 88000];
Mm=[2 2 2 1 2 1 2 2 2 2 2 2 2];
Mk=[3 3 3 2 4 4 4 4 4 4 4 4 4];
Mdelt=[0.1 0.1 0.1 0.1 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05];

 Msavedir={ '/media/09004e3f-3e8e-409e-b718-5a0ca5495abd/home/cenna/Results/2D/Fric800s2_m1_k4/0/' };
Mstart=[ 2000];
Msteps=[2000 ];
Mend=[98000];
Mm=[1 ];
Mk=[4];
Mdelt=[0.01];
figure
for simi=[1] 
savedir=Msavedir{simi};
delt=Mdelt(simi);
delx=1;
start=Mstart(simi);
steps=Msteps(simi);
ends=Mend(simi);
L=1;
m=Mm(simi);
kappa=Mk(simi);
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
Mtn=[start:steps:ends];
try
for tni=1:size(Mtn,2)
    tn=Mtn(tni);
    grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    grains(grains<10)=[];
    R=(grains/pi).^0.5;
    Rbar(tni)=mean(R)*delx;
    numgrains(tni)=length(grains);
    biggrain(tni)=grains(1)*delx^2;
end
end
Mtn=[start:steps:start+steps*(length(Rbar)-1)];
timevec=Mtn*delt;
timevec=timevec-timevec(1);
tau=timevec*intenergy*mobility;
wwwww=intenergy*mobility

% figure
plot(timevec,Rbar,'o')
ylabel('Mean Radius of Grains');xlabel('Time')
%  figure
 hold on
% plot(tau,Rbar.^2,'bd','LineWidth',1.5)
% ylabel('<R>^2','FontSize',14);xlabel('\tau','FontSize',14);grid on;
% hold on

pp=polyfit(tau, Rbar.^2,1);
slope(simi)=pp(1);
end
% set(gca,'xlim',[0 6000])
set(gca,'Box', 'on')
set(gca,'FontSize',12,'LineWidth',2)
mean(slope)
%% beta from kinetics
timevec(80:end)=[];
Rbar(80:end)=[];
pp=polyfit(timevec,Rbar.^2,1);
slope_from_the_graph=pp(1)
L=1;
m=1;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
beta=slope_from_the_graph/2/mobility/intenergy
%% beta from rates
L=1;
m=2;
kappa=3;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
spline1 = spaps(timevec,Rbar,10);
SpRbar=fnval(spline1,timevec);
% plot(timevec,Rbar,'.',timevec,SpRbar,'r-')
dRdt=fnval(fnder(spline1,1),timevec);
beta=dRdt.*SpRbar/mobility/intenergy;
% figure
plot(Rbar,beta,'o')
%% Growth rate vs. size

figure
mkdir([savedir 'Rate_RA/'])
tni=20
for tni=1:size(GrainStat,2)-1
    tn=Mtn(tni);
    grains=GrainStat{tni};
    grains1=GrainStat{tni+1};
    nonzerogi=find(grains1>10);
    grainsvol=grains(nonzerogi)*delx^2;
    grainsvol1=grains1(nonzerogi)*delx^2;
    RA=(1/pi)^(1/2)*grainsvol.^(1/2);
    RA1=(1/pi)^(1/2)*grainsvol1.^(1/2);
    GrowthRate=(RA1-RA)/(steps*delt);
    plot(RA/mean(RA),GrowthRate.*RA,'.')
%     axis([0 2 -1 1])
    title(['Time= ' num2str(tn)])
    ylabel('Growth Rate (dR_A/dt) \times Grain Radius (R_A)','FontSize',14)
    xlabel('Grain Equivalent Radius','FontSize',14)
    title(['Time Step= ' num2str(tn)],'FontSize',14);grid on
    set(get(gcf,'Children'),'LineWidth',1.2)
    set(get(gcf,'Children'),'FontSize',14)
    print([savedir 'Rate_RA/' num2str(tni) '.png'],'-dpng','-r200',gcf)
end

%% Growth rate average
clear
savedir='/media/Disk2/Results/2D/Fric2000_Pz0_m1_k2_run2/'
delt=0.1;
delx=1;
start=2000;
steps=50;
L=1;
m=4;
kappa=8;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
Mtn=[start:steps:30000];
for tni=1:size(Mtn,2)-1
    tn=Mtn(tni);
    tn1=Mtn(tni+1);
    grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    grains1=importdata([savedir 'GrainStat_' num2str(tn1) '.txt']);
    nonzerogi=find(grains1>10);
    grainsvol=grains(nonzerogi)*delx^2;
    grainsvol1=grains1(nonzerogi)*delx^2;
    RA=(1/pi)^(1/2)*grainsvol.^(1/2);
    RA1=(1/pi)^(1/2)*grainsvol1.^(1/2);
    GrowthRate=(RA1-RA)/(steps*delt);
    GrowthRateBar(tni)=mean(GrowthRate.*RA);
    Rbar(tni)=mean(RA);
end
   
beta=0.2
dRbardt=beta*mobility*intenergy./Rbar;
plot(dRbardt,GrowthRateBar,'.','MarkerSize',8)
xlabel('Growth Rate (d<R>/dt)','FontSize',14)
ylabel('<dR/dt>','FontSize',14)
set(gca,'LineWidth',2)
set(gca,'FontSize',12)


%% size histograms from grain stats
clear
 savedir='/home/cenna/Results/2D/Fric800s2_m1_k4/0';
 tn=20000
Rbar=HistGrainSize(savedir,tn)
 savedir='/home/cenna/Results/2D/Fric1000s2_m2_k3_phi_init1000/20_2/';
 tn=92000
Rbar=HistGrainSize(savedir,tn)
 %% Curvature Maps
 savedir='/home/cenna/Results/2D/Fric800s2_m1_k4/0/';
 tn=20000
 mboxsize=1600;
 Khist=CurvatureSurf(savedir,tn,mboxsize);
%% Curvature Histograms
savedir='/home/cenna/Results/2D/Fric800s2_m1_k4/0/';
tn=2000
m=1
kappa=4
Pz=0
maxK=0.025;
dosave='nosave';
PlotType='plot';

Khist=HistCurvature(savedir,m,kappa,Pz,tn,maxK,dosave,PlotType);
%% Face numbers analysis
savedir='/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz20/';
start=2000
steps=2000
ending=10000
FaceNumberCalculator(savedir,start,steps,ending)

%% number of fasces vs. size
clear
figure
delx=1;
delt=0.05;
savedir=['/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz10_run2/'];
start=2000
steps=2000
Mtn=[start:steps:80000];
for tni=1:size(Mtn,2)
    tn=Mtn(tni)
    load([savedir '/facenums/facenum_' num2str(tn) '.mat'])
    grainsvol=grains(nonzerogi)*delx^2;
    R=(grainsvol/pi).^0.5;
    Rbar=mean(R);
    plot(facenum(1:length(nonzerogi)),R/Rbar,'.','MarkerSize',8)
    ylabel('Normalized Grains Equivalent Radius (R/<R>)','FontSize',14)
    xlabel('Number of Faces (N)','FontSize',14)
    title(['Time Step= ' num2str(tn)],'FontSize',14);grid on
    set(gca,'LineWidth',1.2)
    set(gca,'FontSize',12)
    axis([0 12 0 2.5])
%     print([savedir '/facenums/R_N_' num2str(tni) '.png'],'-dpng','-r200',gcf)
   pause(0.2)
end

%% Von-Nuemann Mullins Relation
clear
figure
delx=1;
delt=0.05;
savedir='/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz10/';
L=1;
m=2;
kappa=4;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
mkdir([ savedir '/facenums/von_mull/'])
start=4000
steps=2000
statstep=100
Mtn=[start:steps:900000];
Mtn=4000
for tni=1:size(Mtn,2)
    tn=Mtn(tni);
    tnp1=Mtn(tni)+statstep;
    tnm1=Mtn(tni)-statstep
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

    GrowthRate=(RAp1-RAm1)/(2*statstep*delt);
    GrowthRateBar(tni)=mean(GrowthRate.*RA);
    Rbar(tni)=mean(RA);
    load([savedir 'facenums/facenum_' num2str(tn) '.mat']);
    facenum=facenum+1;
    VonYClass=[1:12]*nan;
    for FaceClass=3:12
        items=find(facenum==FaceClass);
        GRClass=GrowthRate(items);
        RAClass=RA(items);
        VonYClass(FaceClass)=mean(GRClass.*RAClass);
    end
    plot(facenum,GrowthRate.*RA/mobility/intenergy,'.','MarkerSize',8)
    ylabel('(dR/dt) \times R / M / \sigma','FontSize',14)
    xlabel('Number of Faces (N)','FontSize',14)
    title(['Time Step= ' num2str(tn)],'FontSize',14);grid on
    set(gca,'LineWidth',2)
    set(gca,'FontSize',14)
    axis([3 12 -1.0 1.0])
    vonY=1/6*([2:12]-6);
    hold on
    plot([2:12],vonY,'r','LineWidth',1.5)

    %      pp=polyfit(facenum,GrowthRate.*RA,1);
    FaceClassi=[1:12];
    FaceClassi(isnan(VonYClass))=[];
    VonYClass(isnan(VonYClass))=[];
    pp=polyfit(FaceClassi,VonYClass/mobility/intenergy,1);
    voni=polyval(pp,[2:12]);
%     plot([2:12],voni,'g','LineWidth',1)
    plot(FaceClassi,VonYClass/mobility/intenergy,'ks')
    print([savedir '/facenums/von_mull/von_mull_' num2str(tni) '.png'],'-dpng','-r200',gcf)
    pause(0.5)
    hold off
end



%% Rc vs Pz
% m2 k3
MPz=[50 30 25 20 20 15 10 ]*10^-3;
MRc=[378.3 986.5 1713 2225 2546 3414 7100].^0.5;
%% init 1000
MPz=[50 30 25 20 15 10 ]*10^-3;
MRc=[19.45 26.07 41.39 46.7 58.81 78.38];
%% init 5000
MPz5000=[30 20 15]*10^-3;
MRc5000=[26.15 31.44 45.66];


% m1 k2
MPz=[80 50 30 20]*10^-3;
MRc=[757 984 1714 2600].^0.5;

figure
plot(1./MPz,MRc,'s');grid on
xlabel('1/Pz')
ylabel('R_{lim}')
%%
 figure
 plot(log(MPz),log(MRc),'rs')
xlabel('Ln(Pinning Pressure), Pz')
ylabel('Ln(Limiting Grain Size), Ln(R_c)');grid on;box on

%% calulating alpha from limiting grain size
alpha=MPz./intenergy.*MRc;

figure
plot(MPz, alpha,'d');grid on;box on
xlabel('Pinning Pressure, Pz')
ylabel('Alpha')

%% analytical comparison
L=1;
m=2;
kappa=4;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
% slope of Equvalent r vs time is:
slope=2*mobility*1*intenergy
slope_from_the_graph=0.7;
alpha=slope_from_the_graph/slope


Pz=0.02;
Rc=2*alpha*intenergy/Pz
Rc^2

alpha=0.848/intenergy
Pzi=linspace(0.005,0.15,100);
Rci=alpha*intenergy./Pzi;
hold on
plot(1./Pzi,Rci)

plot(Pzi,Rci)

%% realistic values for pinning
ezplot('2/3*0.015*r',[0.01 0.1])
hold on
ezplot('2/3*0.045*r',[0.01 0.1])
grid on
ylabel('Particles Volume Fraction (f_v)','FontSize',14)
xlabel('Particle Radius (r)','FontSize',14)
grid on
set(get(gcf,'Children'),'LineWidth',1.2)
set(get(gcf,'Children'),'FontSize',14)
%%
figure
ezcontour('(1/(4/3*pi*r^3))*f', [0 3e-3], [0.01 0.1], 50)
shading interp
xlabel('Particles Volume Fraction (f_v)','FontSize',14)
ylabel('Particle Radius (r)','FontSize',14)
title('Number Density of Particles (N_v)' ,'FontSize',14);grid on
set(get(gcf,'Children'),'LineWidth',1.2)
set(get(gcf,'Children'),'FontSize',14)

% view([0 0])
%%-------------------------------------------------------------------------
%% Grain Statistic with For loop extracted from Inds
clear
mboxsize=1000;
nboxsize=mboxsize;

nuclein=abs(mboxsize*nboxsize/200);
delx=2;
delt=0.15;
steps=500;
savedir='/media/disk/sim_res/2D/Fric1000_m1_k2_2/50/';
Mtn=[steps:steps:200000];
try
    for tni=1:size(Mtn,2)
        tn=Mtn(tni)
        Inds=importdata([savedir 'Inds_' num2str(tn) '.txt']);
        % obtaining statistics
        V=zeros(1,nuclein+1);
        for ni=[1:nuclein+1]
            % search for ni in the Inds matrix
            indsni=find(Inds==ni);
            V(ni)=length(indsni);
        end
        %     V(V<1^3)=[];
        Vbar=mean(V);
        MV{tni}=V;
        MVbar(tni)=Vbar;
    end
catch
end
timevec=Mtn*delt;
plot(timevec,MVbar,'o')
figure
plot(timevec,MVbar.^(2/3),'x')



%% Grain Size from Phi Fields
clear
dirstring='/media/disk/sim_res/2D/Normal_m1_k2/';
ni=0;
delx=2;
delt=0.1;
steps=500
Mtn=[steps:steps:200000];
try
for tni=1:size(Mtn,2)
    tn=Mtn(tni)
    phi=importdata([dirstring 'Fullres_' num2str(tn) '.txt']);
    [areaG,diamG,perimG]=grainstat_simple(phi);
    A{tni}=areaG*delx^2;
    D{tni}=diamG*delx;
    P{tni}=perimG*delx;
    %     M(i)=getframe;
end
catch;end
save(strcat(dirstring,'/','Grains_phi','.mat'))
%% Total grains area
% clear
clear
dirstring='/home/cenna/Results/2D/Fric1000_m1_k2_2/80/';
load(strcat(dirstring,'Grains_phi','.mat'))
delx=2;
delt=0.15;
steps=500;
Mtn=[steps:steps:200000];
try
for tni=1:size(Mtn,2)
        Ai=A{tni};Di=D{tni};Pi=P{tni};
        % removing small dots
        [removeindex]=find(Ai<5);
        Ai(removeindex)=[];
        Di(removeindex)=[];
        Pi(removeindex)=[];
        An(tni)=mean(Ai);
        Dn(tni)=mean(Di);
        Pn(tni)=mean(Pi);
        GBArea(tni)=sum(Pi)/2;
end
catch;end
Mtn=[steps:steps:tn-steps];
timevec=Mtn*delt;
SA=(GBArea/(size(phi,1)*size(phi,2)*delx^2));
%%
figure
hold on
plot(timevec,An/pi)
ylabel('(Equivalent Radius)^{2}');xlabel('Time');grid on
%%
figure
hold on
plot(timevec,1./SA.^2)
ylabel('Inverse Total Grain Boundary Area');xlabel('Time');grid on
%%
% interface area/ grain size, we want to see if alpha is different when we
% add pinning and this ratio may give an insight. Also its good for
% checking self similarity.
figure
plot(timevec, (An/pi).*SA.^2)
ylabel('Total Grain Boundary Area*Grain Size');xlabel('Time');grid on
%%
figure
hold on
plot(timevec, (Dn).*SA)
ylabel('Total Grain Boundary Area*Grain Size');xlabel('Time');grid on
%%
% figure
hold on
plot(1./Dn,SA)
ylabel('Total Grain Boundary Lines in Unit Area (1/m)');xlabel('Inverse Average Grain Diameter (1/m)');grid on
%%
MPz=[80 50 30 20 10 5]
MIAt=[7.249e-6 8.266e-6 1.114e-5 1.4e-5 1.731e-5 2e-5]
figure
loglog(MPz,MIAt,'s')
xlabel('Pinning Pressure, Pz')
ylabel('Inverse Total Area, 1/A_t')
pp=polyfit(log(MPz),log(MIAt),1)

%% B vs n plot for solute drag
clear
b= [200 300 400 500 600 1000 1500 2000 2500 3000 3500 4000 4500 5000 10000 20000 50000 100000 200000 400000]
b=400000
figure
for simnum=1:length(b)
    for runnum=[1]
    savedir=['/media/Disk2/Results/2D/solute_m2_k4/Fric2000_a1_b' num2str(b(simnum)) '_m2_k4_run' num2str(runnum) '/']
    delt=0.05;
    delx=1;
    start=2000;
    steps=2000;
    Mtn=[start:steps:150000];Rbar=[];
    try
        for tni=1:size(Mtn,2)
            tn=Mtn(tni);
            GrainStat{tni}=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
            grains=GrainStat{tni};
            grains(grains<5)=[];
            R=(grains/pi).^0.5;
            Rbar(tni)=mean(R)*delx;
        end
    catch
    end
    timevec=Mtn(1:length(Rbar))*delt;
    % figure
    plot(timevec,Rbar,'.')
    ylabel('Equivalent Raidus ^2');xlabel('Time')
    hold on
    [n(runnum) ci(:,runnum)]=curve_fit(timevec,Rbar);
    end
    nt(simnum)=mean(n);
    maxn(simnum)=max(n);
    minn(simnum)=min(n);
end

%%
figure
semilogx(b,nt,'o')
hold on
% semilogx(b,ci(1,:),'r:')
% semilogx(b,ci(2,:),'g:')
errorbar(b,nt,nt-minn,nt-maxn,'ro')
% errorbar(b,nt,nt-ci(1,:),ci(2,:)-nt,'ro')
xlabel('Solute Drag Parameter (b)')
ylabel('Grain Growth Time Exponent (n)')


%% ------------------------- 2D CLOSE Friction Stuff ----------------------------

%% Particle- Interface
clear
ni=0
figure
tn=357000;
% for tn=[1000:1000:600000]
ni=ni+1;
savedir='/home/cenna/Results/Np/Fv-1_r20_DelG0.01/';
phi=importdata([savedir 'Fullres_' num2str(tn) '.txt']);
%eta(:,:,1)=importdata([savedir 'Eta_0_' num2str(tn) '.txt']);
%eta(:,:,2)=importdata([savedir 'Eta_1_' num2str(tn) '.txt']);
ppf=double(phi<0.0);

% structure
%  figure
% imshow(phi);
surf(phi);
view ([0 90])
shading flat
axis off
box off
%caxis([0 1]);
colormap jet
shading flat
title(['Time Step= ' num2str(tn)])
set(gca,'DataAspectratio',[1 1 1])
pause(0.01);
% end
% plot(Etn,'.-');grid
%% energy
mboxsize=size(phi,1);
nboxsize=size(phi,2);
delx=0.5;
settings.m=2;
settings.kappa=4;
settings.epsilon=5;
settings.accuracy='low';

ppf=double(phi<0.1);
[ME,E]=calculateE(eta,ppf,mboxsize,nboxsize,delx,settings);

surf(ME);
view ([0 90])
shading flat
axis off
box off
axis equal
colorbar
% set(gca,'DataAspectratio',[1 1 1])
%% analysing volume logs
clear
figure; hold on
savedir='/home/cenna/Results/Np/Fv-1_r20/';
delx=1;
nboxsize=100;
mboxsize=200
DelG=0.10
fv=0.0843


voldata=importdata([savedir 'Vollog' '.log']);
 time=voldata(:,1);
 vol=voldata(:,2);
 Etn=voldata(:,3);
 plot(time,vol)
 figure
 plot(vol,Etn)
    %%
    length=vol/(nboxsize*delx*(1-fv));
     plot(time,length)
    pp=polyfit(time,length,1);
    Mvel=pp(1);
    
xlabel('Time')
ylabel('Volume of upper grain')

mobility=Mvel/DelG
L=1;
m=1;
kappa=2;
mobility_theory=3/2*L*sqrt(2*kappa/m)
intenergy=1/3*sqrt(2*m*kappa);
Pz=DelG-Mvel/mobility_theory
%%
figure
Mfv=[0.128 0.113 0.089 0.046];
MPz=[0.0184 0.0157 0.0140 0.0068];
plot(Mfv,MPz,'o')
hold on
fvi=linspace(0.03,0.15,10);
Pzi=3/2*intenergy*fvi/8;
plot(fvi,Pzi,'r')
%%
figure

Mcurvature=1./(Mr*delx);
plot(Mcurvature,Mvel,'o')
xlabel('Intreface Curvature');
ylabel('Interface Velocity');
grid on

% mobility:
sigma=1;
mobility=polyfit(Mcurvature*sigma,Mvel,1)

%% analysing volume logs for CIRCLE grain

clear
figure; hold on
savedir='/home/cenna/Results/2Dclose/diamond_Fric/Pz0/';

Mr=[100];
steps=[2];
tol=0.00000001
L=1;
m=4;
kappa=8;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
alpha=1;
for n=1:size(Mr,2)
    voldata=importdata([savedir 'vollog_50' '.log']);
%     voldata=importdata([savedir 'vollog' '.log']);
    time=voldata(:,1);
    vol=voldata(:,2);
    radius=sqrt(vol/pi);
%     figure
%     plot(time,radius); xlabel('Time'); ylabel('Radius')
    spline1 = spaps(time(1:steps(n):end),radius(1:steps(n):end),tol(n),3);
    Mvel=fnval(fnder(spline1,1),time(1:steps(n):end));
% xlabel('Time')
% ylabel('radius of circular grain')
% figure
Mvel=-Mvel;
Mcurvature=alpha./radius(1:steps(n):end);
Mcurvaturei=linspace(min(Mcurvature), max(Mcurvature),100);
Mveli=interp1(Mcurvature,Mvel/mobility/intenergy,Mcurvaturei);
plot(Mcurvaturei,Mveli,'.')
xlabel('P_d/ \sigma_{gb}');
ylabel('V / (M \sigma_{gb})');
grid on
end
 axis([0 0.18 0 0.18]);box on
hold on
plot([0 0.5],[0 0.5],'r')
title(savedir)

%% Analytical V-DelG for Solute Drag
L=1;
m=2;
kappa=4;
a=1;
b=20000;
syms v;
Pf=a*v/(1+b*v^2);
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
delG=linspace(0,0.05,100);
ni=0
for dG=delG
    ni=ni+1;
    e1=mobility*(dG-Pf)-v;
    e1sol=solve(e1);
    vel_sol(ni)=eval(e1sol(1));
end
figure
hold on
% plot(delG,vel_sol)
% plot([0 max(delG)], [0 max(delG)*mobility],'r')
plot(delG/intenergy,vel_sol/mobility/intenergy,'k')
xlabel('\Delta G/ \sigma_{gb}');
ylabel('V / (M \sigma_{gb})');
plot([0 5],[0 5],'r')
axis([0 max(delG)/intenergy 0 max(delG)/intenergy]);box on
%%
figure
Pfi=subs(Pf,'v',vel_sol);
plot(vel_sol,Pfi)



%% 
figure
plot(time(1:steps(n):end), Mcurvature,'.')
ylabel('\Delta G/ \sigma_{gb}')
xlabel('time')
%% plot quality of time-radius spline fit
figure
plot(time,radius)
xlabel('time')
ylabel('Radius')
hold on
spline1 = spaps(time(1:steps(n):end),radius(1:steps(n):end),tol(n),3);
plot(time(1:steps(n):end),fnval(spline1,time(1:steps(n):end)),'r-')
grid on


%% analysing volume logs for DOME grain

clear
figure; hold on
savedir='/home/cenna/Results/2Dclose/dome_solute/dome_a0.4_b4_m2_k2_s3/';

Mr=[5 10 15 20 30];
delx=2/3;
L=1;
m=2;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
alpha=1.5;
for n=1:size(Mr,2)
%     voldata=importdata([savedir 'vollog_' num2str(r) '.log']);
    voldata=importdata([savedir '/' num2str(Mr(n)) '/' 'vollog' '.log']);
    time=voldata(:,1);
    vol=voldata(:,2);
    domel=(vol-pi*Mr(n)^2/2)/(2*Mr(n));
%     plot(time,domel)
%     xlabel('Time')
%     ylabel('radius of circular grain')
%     figure
Mveln=polyfit(time,domel,1);
Mvel(n)=-Mveln(1);
end

Mcurvature=alpha./Mr;
plot(Mcurvature,Mvel/mobility/intenergy,'o')
xlabel('\Delta G/ \sigma_{gb}');
ylabel('V / (M \sigma_{gb})');
grid on
axis([0 0.5 0 0.5]);box on
hold on
plot([0 0.5],[0 0.5],'r')
title(savedir)

%% volume logs for ONSET of circle grain simulations

clear
figure; hold on
savedir='/home/cenna/Results/2Dclose/diamond_Fric/Onset_m2_k2/';
delx=1;
scale=2;
Mr=[ 100 50 45 40 35];
Mr=50
steps=[100 50 20 1 1 1];
tol=0.01*[1 1 1 1 1];
MPz=[];
L=1;
m=2;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
alpha=0.321;
for n=1:size(Mr,2)
    r=Mr(n);
    voldata=importdata([savedir 'vollog_' num2str(r) '.log']);
    Pz=voldata(:,3);
    movinginds=find((Pz==Pz(end))==1);
    time=voldata(movinginds,1);
    vol=voldata(movinginds,2);
    radius=sqrt(vol/pi);
   % radiusi=linspace(min(radius),max(radius),1000);
   % timei=interp1(time,radius,radiusi);
   % plot(time,radius)
    spline1 = spaps(time(1:steps(n):end),radius(1:steps(n):end),tol(n),3);
    Mvel=fnval(fnder(spline1,1),time(1:steps(n):end));
% xlabel('Time')
% ylabel('radius of circular grain')

Mvel=-Mvel;
Mcurvature=alpha./radius(1:steps(n):end);
inds=find(Mcurvature>0.035);
Mcurvature(inds)=[];
Mvel(inds)=[];
plot(Mcurvature,Mvel/mobility/intenergy)
xlabel('\Delta G/ \sigma_{gb}');
ylabel('V / (M \sigma_{gb})');
grid on
MPz=[MPz Pz(end)];
end
Normalized_Pz=MPz/intenergy
axis equal
axis([0 0.035 0 0.035]);box on; 
plot([0 0.035], [0 0.035],'r')
%%
observedsrart=[0.004570 0.01151 0.02247]
Normalized_Pz=[0.0040    0.0104    0.0206]
figure
plot(Normalized_Pz,(observedsrart-Normalized_Pz)./Normalized_Pz,'s-')
grid on
xlabel('Normalized Driving Pressure')
ylabel('(observedsrart-Normalized Pz)/Normalized Pz')

figure
plot(Normalized_Pz,observedsrart,'s-')
grid on;box on; hold on
plot([0 max(Normalized_Pz)] , [0 max(Normalized_Pz)],'r')
xlabel('Normalized Driving Pressure')
ylabel('Observed Starting Pressure')
title(savedir)

%% analytical exprression for solute drag
syms v r M a b sigma;
Pf=a*v/(1+b*v^2);
e=v-M*(sigma/r-a*v/(1+b*v^2))
velocity=solve(e,'v')
velocity=simple(velocity)
f=velocity(1);

L=1;
m=4;
kappa=8;
a=1;
b=1000;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
M=mobility;
sigma=intenergy;
f=eval(f)
t=int(1/f,'r')

%% ----------------------------------------- Particles N   -----------------------


%% particlesN 2D
clear
savedir='/home/cenna/Results/Np/';
delx=1;
mboxsize=100;
DelG=0.05;

voldata=importdata([savedir 'Vollog' '.log']);
fv=voldata(1,1);
r=voldata(1,2);
time=voldata(2:end,1);
vol=voldata(2:end,2);
length=vol/(mboxsize*delx*(1-fv));
pp=polyfit(time,length,1);
Mvel=pp(1);

figure; hold on
plot(time,length)
xlabel('Time')
ylabel('Volume of upper grain')

% mobility:
Meff=Mvel/DelG;
L=1;
m=4;

kappa=8;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
Pftot=DelG-Mvel/mobility
%% Pz analytic
Pz_analitic=3/2*intenergy*fv/r
%% analytical exprression for solute drag of flat interface
L=1;
m=1;
kappa=2;
a=1;
b=1000;

mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);

vi=linspace(0,1,1000);
Ps=a*vi./(1+b*vi.^2);

hold on
plot(vi,Ps)
xlabel('V');ylabel('P_s')
%%
syms v ;
e=v-mobility*(G-a*v/(1+b*v^2));
velocity=solve(e,'v');
velocity=eval(velocity)
velocity=Mvel
Ps=a*velocity(1)/(1+b*velocity(1)^2)
Pz=0.0379;
error=(Pftot-(Ps+Pz))/Pftot*100

%% ----------------------------------------- Particles N Solute  -----------------------
%% analysing volume logs
clear
L=1;
m=2;
kappa=4;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);

MG=[ 350 300 250 200 150 100 80 50 40 30 20 10]
MDelG=MG/1000;
for ni=1:length(MG)
    savedir=['/home/magnetadmin/Documents/Results/Np_Solute/Np1000_m2_k4_r8/' num2str(MG(ni)) '/'];
delx=0.5;
nboxsize=1000;
DelG=MG(ni)/1000;
voldata=importdata([savedir 'Vollog' '.log']);
fv=voldata(1,1);
    time=voldata(2:end,1);
    vol=voldata(2:end,2);
    l=vol/(nboxsize*delx*(1-fv));
    pp=polyfit(time,l,1);
    Mvel(ni)=pp(1);
    %figure
% plot(time,l)
% xlabel('Time')
% ylabel('Volume of upper grain')
end
% figure; hold on
% plot(MDelG,Mvel,'d','LineWidth',1.5); grid on
% xlabel('\Delta G','FontSize',14);ylabel('v','FontSize',14)
% set(gca,'LineWidth',2,'FontSize',13)

figure
MPf=MDelG-Mvel/mobility
plot(Mvel,MPf,'s')
xlabel('V','FontSize',14);ylabel('P_f','FontSize',14)
set(gca,'LineWidth',2,'FontSize',13)
%%
mobility=Mvel/DelG
L=1;
m=1;
kappa=2;
mobility_theory=3/2*L*sqrt(2*kappa/m)
intenergy=1/3*sqrt(2*m*kappa);
Pz=DelG-Mvel/mobility_theory


