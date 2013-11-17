%% Visualization
%blah blah

%% 3D slice full reading (phi)
clear
tn=200
ni=20;
figure
% for tn=[20:20:64000]
ni=ni+1;
savedir='/home/cenna/Results/Rec/';
phidata=importdata([savedir 'Fullres_' num2str(tn) '.txt']);
mboxsize=size(phidata,2);nboxsize=sqrt(size(phidata,1));lboxsize=nboxsize;
% x=[0:nboxsize-1]*delx;/home/cenna/Results/test/EngDen_900.txt
% y=[0:mboxsize-1]*delx;
% z=[0:lboxsize-1]*delx;
% [x,y,z]=meshgrid(x,y,z);

phi=zeros(mboxsize,nboxsize,lboxsize);
for m=1:mboxsize
    for n=1:nboxsize
        for l=1:lboxsize
            phi(m,n,l)=phidata(m+(l-1)*mboxsize,n);
        end
    end
end
clear phidata;

% figure; 
   slice(phi,[1 (mboxsize)],[1 (nboxsize)],[1 (lboxsize)]);colormap gray; shading flat; axis equal; box on; title (['time step = ' num2str(tn) ]) % colorbar; title ([savedir ' (tn=' num2str(tn) ')'])
% slice(phi,[ (mboxsize)],[ (nboxsize)],[1 (lboxsize)/2]);colormap gray; shading flat; axis equal; box on; title (['time step = ' num2str(tn) ]) % colorbar; title ([savedir ' (tn=' num2str(tn) ')'])
 % axis([1 (mboxsize) 1 (nboxsize) 1 (lboxsize)])
%  surf(phi(:,:,(lboxsize)/2)); view([0 90]);colormap gray; shading flat; axis equal; box on; title (['time step = ' num2str(tn) ])
 % mkdir([savedir 'phi/']);
%  print([savedir 'phi/' num2str(ni) '.png'],'-dpng','-r200',gcf)
%   end
sliceomatic(phi)
% plot([0.5:0.5:35], phi(:,1,1))

%% Fourer Transform
Y=fftn(phi); % phi is a grayscale image in 3D array. 
Y=(real(Y).^2+imag(Y).^2);
Y(1,1,:)=0;
Y(:,1,1)=0;
Y(1,:,1)=0;
figure
cut=30 % frequencies larger than 30 [1/pixcel] is removed because nothing is there.
slice((Y(1:cut,1:cut,1:cut)),[1 (cut)],[1 (cut)],[1 (cut)]);
colormap jet
shading interp
axis equal
box on
colorbar
xlabel('kx');ylabel('ky');zlabel('kz');

 
%% 3D slice full reading (INDS)
clear
tn=120000;
ni=1;
savedir='/media/Disk2/Results/3D/Fric300_Pz0.010_m1_k2_init4000_run1/';
% figure
% for tn=[100:100:1900]
phidata=importdata([savedir 'Inds_' num2str(tn) '.txt']);
mboxsize=size(phidata,2);nboxsize=mboxsize;lboxsize=mboxsize;
Inds=zeros(mboxsize,nboxsize,lboxsize);
for m=1:mboxsize
    for n=1:nboxsize
        for l=1:lboxsize
            Inds(m,n,l)=phidata(m+(l-1)*mboxsize,n);
        end
    end
end
clear phidata;
figure;
slice(Inds,[1 (mboxsize)],[1 (nboxsize)],[1 (lboxsize)]);
axis([1 (mboxsize) 1 nboxsize 1 lboxsize]);
colormap jet;
shading flat
axis equal;
box on;
colorbar;
title ([savedir ' (tn=' num2str(tn) ')'])
 mkdir([savedir 'Inds/' ])
 print([savedir 'Inds/' num2str(ni) '.png'],'-dpng','-r200',gcf)
ni=ni+1;
%  end
%sliceomatic(Inds)
%% finding position of a particular grain

grainind=7;
gi=find(Inds==grainind);
[a,b,c]=ind2sub(size(Inds),gi);
if ~isempty(a)
    figure
    plot3(a,b,c,'.')
    axis([0 size(Inds,1) 0 size(Inds,2) 0 size(Inds,3)]);grid on; box on
end

%% Grain Size from GrainStat files
clear
savedir='/home/cenna/Results/Rec/';
delt=0.05;
delx=1;
start=100; steps=100;
Mtn=[start:steps:2000000];
try
    for tni=1:size(Mtn,2)
        tn=Mtn(tni);
        GrainStat=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
        grains=GrainStat;
        [a,b]=max(GrainStat);
        abnormG(tni)=((3/4/pi)^(1/3))*grains(b).^(1/3);
        grains(grains<3^3)=[];
        radius=((3/4/pi)^(1/3))*grains.^(1/3);
        Rbar(tni)=mean(radius);
        numgrais(tni)=length(grains);
        bigR(tni)=radius(1);
        bigV(tni)=grains(1);
    end
catch
end     
Mtn=[start:steps:tn-steps];
timevec=Mtn*delt;
%timevec=timevec-timevec(1);
%figure
hold on
plot(timevec,Rbar,'r.')
ylabel('Equivalent Radius');xlabel('Time');grid on
% 
% figure
% plot(timevec,Rbar/Rbar(end))
%% beta
L=1;
m=1;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
beta=diff(Rbar).*Rbar(1:end-1)/(steps*delt)/mobility/intenergy;
figure
plot(Rbar(1:10),beta(1:10),'o')
%% Grain Size Distribution
clear
savedir='/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.015_m1_k2_init2000/';
delt=0.1;
delx=2;
start=2000; steps=2000;
ending=8000000;
Mtn=[start:steps:ending];
ni=0;
figure
for tni=1:size(Mtn,2)
    tn=Mtn(tni);
    GrainStat{tni}=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    grains=GrainStat{tni};
    grains(grains<25)=[];
    grains=grains*delx^3;
    radius=((3/4/pi)^(1/3))*grains.^(1/3);
    Rbar=mean(radius);
    % [a,b]=hist(log(grains)/log(mean(grains)),40);
    ri=linspace((2.5), (5.5),20);
    [a,b]=hist(log(radius),ri)
    bar(b,a/length(radius));
    hold on;
    plot([log(Rbar) log(Rbar)],[0 100],'r','LineWidth',2);
%     plot([log(Rbar/Rbar) log(Rbar/Rbar)],[0 100],'r','LineWidth',2);  
    xlabel('log(Grain Equivalent Radius)','FontSize',14)
    ylabel('Counts','FontSize',14)
    hold off;
    set(gca,'LineWidth',1.5)
    set(gca,'FontSize',14)
    axis([min(ri) max(ri) 0 0.25])
%     axis([-2 1 0 0.25])
    title (['time =' num2str(tn*delt) ])
    ni=ni+1;
%     mkdir([savedir 'Ghist/']);
%     print([savedir 'Ghist/' num2str(ni) '.png'],'-dpng','-r200',gcf)
%     clf
 pause(0.1)
end


%% Rc vs Pz
MPz=[0.01 0.015 0.018 0.020 0.025 0.03]
MRc=[41.53 29.95 22.84 20.8 18.73 13.38]
 figure
 plot(1./MPz,MRc,'ro')
xlabel('1/Pz')
ylabel('R_{lim}');grid on;box on

%% analytical comparison
L=1;
m=1;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
% slope of Equvalent r vs time is:
slope=2*mobility*1*intenergy
slope_from_the_graph=0.9988
alpha=slope_from_the_graph/slope

%% Rc vs Pz analytical
m=2;
kappa=3;
intenergy=1/3*sqrt(2*m*kappa);

alpha=0.6368;
Pzi=linspace(0.008,0.05,100);
Rci=alpha*intenergy./Pzi;
hold on
plot(1./Pzi,Rci)
% plot(Pzi,Rci)

%% Rc vs Pz Simulation
Pz=[0.05 0.04 0.03];
Rc=[195 657 3200 ].^0.5;
figure
loglog(Pz,Rc,'s')
polyfit(log(Pz),log(Rc),1)
%%-------------------------------------------------------------------------


%% Finding neighbors of a grain (number of faces)
clear
mboxsize=512;
nboxsize=mboxsize;
lboxsize=mboxsize;
nuclein=abs(mboxsize*nboxsize*lboxsize/1000);
delx=1;
delt=0.1;
savedir='/home/magnetadmin/Documents/Results/3D/Fric512_Pz0_m1_k2_init0/';
steps=1000
start=1000
Mtn=[start:steps:200000];
% Mtn=[50000];
tn=8000;
% finding the bigest grain index:
grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
[a,ind]=max(grains);grains(ind)=0;
[a,ind2]=max(grains);% second largest grain at the end
%%
se=strel(ones(2,2,2));
try
    for tni=1:size(Mtn,2)
        tn=Mtn(tni)
        Inds=importdata([savedir 'Inds_' num2str(tn) '.txt']);
        BWInds=uint8(zeros(size(Inds)));
        % obtaining statistics
        nii=1;
        for ni=ind
            % search for ni in the Inds matrix
            BWInds(Inds==ni)=1;
            BWInds=imdilate(BWInds,se);
            Inds2BW=Inds.*double(BWInds);
            listi=Inds(Inds2BW~=0);
            facenum(tni)=0;
            for i=1:length(listi)
                indi=listi(i);
                if indi~=0
                    facenum(tni)=facenum(tni)+1;
                    listi(listi==indi)=0;
                end    
            end
            nii=nii+1;
        end
    end
catch;end

% time history, number of naces vs. time for a grain with index ind
Mtn=[steps:steps:tn-steps];
timevec=Mtn*delt;
figure
plot(timevec,facenum,'.')
title(['Number of faces for Grain with index= ' num2str(ind)])

%% Topology Distributions
clear
Msavedir={['/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.010_m1_k2_init4000_run1/']};
Mmboxsize=[300];
Mstart=[4000] ;
Msteps=[2000];
Mends=[80000];
for simi=1
    mboxsize=Mmboxsize(simi);
    nboxsize=mboxsize;
    lboxsize=mboxsize;

    delx=1;
    delt=0.1;
    savedir=Msavedir{simi};
    mkdir([savedir '/facenums'])
    steps=Msteps(simi);
    start=Mstart(simi);
    ends=Mends(simi);
    Mtn=[start:steps:ends];
    %  Mtn=[1000];
    se=strel(ones(2,2,2));
    for tni=1:size(Mtn,2)
        tn=Mtn(tni)
        Inds=importdata([savedir 'Inds_' num2str(tn) '.txt']);
        % finding index of grains with non zero size:
        grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
        [nonzerogi]=find(grains>10);
        % obtaining statistics
        facenum=zeros(size(nonzerogi));
        for gi=1:length(nonzerogi)
            ni=nonzerogi(gi);
            % search for ni in the Inds matrix
            BWInds=zeros(size(Inds),'uint8');
            BWInds(Inds==ni)=1;
            BWInds=imdilate(BWInds,se);
            Inds2BW=Inds.*double(BWInds);
            listi=Inds(Inds2BW~=0);
            facenum(gi)=-1;
            for i=1:length(listi)
                indi=listi(i);
                if indi~=0
                    facenum(gi)=facenum(gi)+1;
                    listi(listi==indi)=0;
                end
            end
        end
        save([savedir '/facenums/facenum_' num2str(tn) '.mat'], 'facenum','grains', 'nonzerogi','tn')
    end
end

%%
figure
delx=2;
delt=0.1;
savedir=['/home/cenna/Results/3D/NormalGG150_m1/19points/'];
steps=200
start=200
Mtn=[start:steps:200000];
for tni=1:size(Mtn,2)
    tn=Mtn(tni)
    load([savedir '/facenums/facenum_' num2str(tn) '.mat'])
    Rm=((3/4/pi)*mean(grains(nonzerogi))*delx^3)^(1/3);
grainsvol=grains(nonzerogi)*delx^3;
    plot(((3/4/pi)^(1/3))*grainsvol.^(1/3)/Rm,facenum(1:length(nonzerogi)),'.','MarkerSize',8)
    xlabel('Normalized Grains Equivalent Radius (R_A/R_m)','FontSize',14)
    ylabel('Number of Faces (N_A)','FontSize',14)
    title(['Time Step= ' num2str(tn)],'FontSize',14);grid on
    set(get(gcf,'Children'),'LineWidth',1.2)
    set(get(gcf,'Children'),'FontSize',14)
    axis([0 2 0 25])
    print([savedir '/facenums/R_NA_' num2str(tni) '.png'],'-dpng','-r200',gcf)

end

%% Mullins Relation
clear
figure
delx=1;
delt=0.1;
L=1;
m=1;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);

savedir='/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.010_m1_k2_init4000_run1/';
% savedir='/home/magnetadmin/Documents/Results/3D/Fric512_Pz0_m1_k2_init0/';
mkdir([ savedir '/facenums/von_mull/'])
start=40000
steps=40000
ends=40000
Mtn=[start:steps:ends];
grainstep=200;
for tni=1:size(Mtn,2)
    tn=Mtn(tni);
    tn1=Mtn(tni)+grainstep;
    load([savedir '/facenums/facenum_' num2str(tn) '.mat']);
    grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    grains1=importdata([savedir 'GrainStat_' num2str(tn1) '.txt']);
    grainsvol=grains(nonzerogi)*delx^3;
    grainsvol1=grains1(nonzerogi)*delx^3;
    RA=((3/4/pi)^(1/3))*grainsvol.^(1/3);
    RA1=((3/4/pi)^(1/3))*grainsvol1.^(1/3);
    GrowthRate=(RA1-RA)/grainstep/delt;
    dVdt=(grainsvol1-grainsvol)/grainstep/delt;
    facenum=facenum+1;
   
%     plot(facenum+1,GrowthRate.*RA,'.','MarkerSize',8)
%     ylabel('Growth Rate (dR_A/dt) \times Grain Radius (R_A)','FontSize',14)
%     xlabel('Number of Faces (N_A)','FontSize',14)
%     title(['Time Step= ' num2str(tn) ', Number of Grains= ' num2str(length(RA))],'FontSize',14);grid on
%     set(gca,'LineWidth',2)
%     set(gca,'FontSize',14)
%     axis([3 25 -5 3])
    
    plot(facenum+1,dVdt.*grainsvol.^(-1/3)/mobility/intenergy,'.','MarkerSize',8)
    ylabel('V^{-1/3}dV/dt/2/M/\sigma','FontSize',14)
    xlabel('Number of Faces (N_A)','FontSize',14)
    title(['Time Step= ' num2str(tn) ', Number of Grains= ' num2str(length(RA))],'FontSize',14);grid on
    set(gca,'LineWidth',2)
    set(gca,'FontSize',14)
   % axis([3 25 -5 3])
    
   % print([savedir '/facenums/von_mull/' num2str(tni) '.png'],'-dpng','-r200',gcf)
%   pause(1)
end
Rbar=mean(RA)
%% analytical Mullins
% Relationship between N and R
figure
plot(facenum,RA,'*')
pp=polyfit(facenum,RA,1)
n=4:35
Ri=polyval(pp,n);
hold on
plot(n,Ri,'c')


%% Hilgenfildt

n=4:35
xn=6*(1-2./n);
Dn=2*atan(sqrt(4*(sin(pi./xn)).^2-1));
GH=6/2^(2/3)*(3/4/pi)^(1/3)*(tan(Dn/2)).^(1/3).*(pi/3-Dn).*((n-2).*tan(pi./xn)).^(2/3);
Y=GH*2;
hold on
plot(n,Y,'k')

%% Hilgenfildt with pinning
L=1;
m=1;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
Pz=0.01;
z=Pz/intenergy;

n=4:35
xn=6*(1-2./n);
Dn=2*atan(sqrt(4*(sin(pi./xn)).^2-1));
GH=2*6/2^(2/3)*(3/4/pi)^(1/3)*(tan(Dn/2)).^(1/3).*(pi/3-Dn).*((n-2).*tan(pi./xn)).^(2/3);

% Rcr=9/8*Rbar
% Ri=Rcr*GH*(4/3*pi)^(-2/3)+Rcr;
pp=polyfit(facenum,RA,1)
% pp=[0.0019   -0.1532    4.8939  -14.2447]
Ri=polyval(pp,n)*4;
    for i=1:length(n)
        Y(i)=0; % else when pinning is stronger
         if (GH(i)+z*Ri(i))<0
             Y(i)=(GH(i)+z*Ri(i));
         end
         if (GH(i)-z*Ri(i))>0
             Y(i)=(GH(i)-z*Ri(i));
         end
    end

hold on
plot(n,Y,'r','LineWidth',1.5)



%% Growth Rate and beta
clear
figure
delx=2;
delt=0.1;
savedir='/media/Disk2/Results/3D/Fric300_Pz0.015_m1_k2_init2000/';
L=1;
m=1;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
% mkdir([ savedir '/Rate/'])
start=2000
steps=100
Mtn=[start:steps:200000];

for tni=1:size(Mtn,2)
    tn=Mtn(tni);
    tn1=Mtn(tni+1);
    grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    grains1=importdata([savedir 'GrainStat_' num2str(tn1) '.txt']);
    nonzerogi=find(grains1>10);
    grainsvol=grains(nonzerogi)*delx^3;
    grainsvol1=grains1(nonzerogi)*delx^3;
    RA=((3/4/pi)^(1/3))*grainsvol.^(1/3);
    RA1=((3/4/pi)^(1/3))*grainsvol1.^(1/3);
    GrowthRate=(RA1-RA)/steps/delt;
    plot(RA/mean(RA),GrowthRate,'.','MarkerSize',8)
    ylabel('Growth Rate (dR_A/dt)','FontSize',14)
    xlabel('Normalized Grain Radius','FontSize',14)
    title(['Time Step= ' num2str(tn)],'FontSize',14);grid on
    set(get(gcf,'Children'),'LineWidth',1.2)
    set(get(gcf,'Children'),'FontSize',14)
    axis([0 3.2 -0.4 0.1])
%     print([savedir '/Rate/Rdot_R' num2str(tni) '.png'],'-dpng','-r200',gcf)
  pause(0.1)
end

%% distributions

figure
tni=40;
grains=(MV{tni});
grains=grains(grains>5^3);
loggrains=log10(grains);
hist(grains/mean(grains))

%% animation
figure
clear
ni=1;
for tn=[3500:1000:120000];
    mboxsize=150;
    nboxsize=mboxsize;
    lboxsize=mboxsize;
    delx=2;
    delt=0.1;
    x=[0:nboxsize-1]*delx;
    y=[0:mboxsize-1]*delx;
    z=[0:lboxsize-1]*delx;
    [x,y,z]=meshgrid(x,y,z);
    savedir='/home/cenna/sim_res/3D/NormalGG150_m1/19points/';
    phidata=importdata([savedir 'Fullres_' num2str(tn) '.txt']);
    phi=zeros(mboxsize,nboxsize,lboxsize);
    for m=1:mboxsize
        for n=1:nboxsize
            for l=1:lboxsize
                phi(m,n,l)=phidata(m+(l-1)*mboxsize,n);
            end
        end
    end
    clear phidata;
    slice(x,y,z,phi,[0 (mboxsize-1)*delx],[0 (nboxsize-1)*delx],[0 (lboxsize-1)*delx]);
    axis([0 (mboxsize)*delx 0 (nboxsize)*delx 0 (lboxsize)*delx])
    colormap gray
    shading flat
    axis equal
    box on
    colorbar
    title (['time =' num2str(tn*delt) ])
    print([savedir 'Inds/' num2str(ni) '.png'],'-dpng','-r200','-f1')
    ni=ni+1;
end

%% animation INDS
figure
clear
ni=1;
savedir='/home/cenna/Results/3D/Fric250_Pz0.01_m1_k2/';
mkdir([savedir 'Inds'])
for tn=[1000:1000:132500];
    mboxsize=250;
    nboxsize=mboxsize;
    lboxsize=mboxsize;
    delx=2;
    delt=0.1;
    x=[0:nboxsize-1]*delx;
    y=[0:mboxsize-1]*delx;
    z=[0:lboxsize-1]*delx;
    [x,y,z]=meshgrid(x,y,z);
    phidata=importdata([savedir 'Inds_' num2str(tn) '.txt']);
    Inds=zeros(mboxsize,nboxsize,lboxsize);
    for m=1:mboxsize
        for n=1:nboxsize
            for l=1:lboxsize
                Inds(m,n,l)=phidata(m+(l-1)*mboxsize,n);
            end
        end
    end
    slice(x,y,z,Inds,[0 (mboxsize-1)*delx],[0 (nboxsize-1)*delx],[0 (lboxsize-1)*delx]);
    axis([0 (mboxsize)*delx 0 (nboxsize)*delx 0 (lboxsize)*delx]);
    colormap jet;
    shading flat
    axis equal;
    box on;
    colorbar;
    title (['time =' num2str(tn*delt) ])
    print([savedir 'Inds/' num2str(ni) '.png'],'-dpng','-r200','-f9')
    ni=ni+1;
end


%% Analytical relationship for kinetics of grain growth with pinning

%% comparison
clear
savedir='/media/Disk2/Results/3D/Fric300_Pz0.013_m1_k2_init4000/';
delt=0.1;
delx=1;
start=4000; steps=100;
Mtn=[start:steps:550000];
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
Mtn=[start:steps:tn-steps];
timevec=Mtn*delt;
figure
hold on
plot(timevec-start*delt,Rbar,'b.')
ylabel('Equivalent Radius');xlabel('Time');grid on

%%
clear
L=1;
m=1;
kappa=2;
Pz=0.013;
alpha= 0.12;
beta=0.20;
R0=19.48;
t0=0;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
alphasigma=alpha*intenergy;
M=mobility;
Rlim=alphasigma/beta/Pz;
syms R t
radiuspart='1/(alphasigma/R-beta*Pz)';
radiuspart=subs(radiuspart,'alphasigma',alphasigma);
radiuspart=subs(radiuspart,'Pz',Pz);
radiuspart=subs(radiuspart,'beta',beta');
Ri=[linspace(R0,Rlim-0.00001,50) linspace(Rlim-0.00001,Rlim-eps,5)];
for ni=1:length(Ri)
    f=int(radiuspart,'R',R0,Ri(ni));
    ti(ni)=eval(f)/M;
end

hold on
plot(ti,Ri)
axis([0 1.5e4 18 32])

%% find alpha and beta to fit the kinetics
clear
R0=19.48;
Rlim=31.70
L=1;
m=1;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
Pz=0.013;
alphai=linspace(0.05,0.6,10);
betai=alphai*intenergy/Rlim/Pz
for curvei=1:length(alphai)
    alpha=alphai(curvei);
    beta=betai(curvei);
    alphasigma=alpha*intenergy;
    M=mobility;
    Rlim=alphasigma/beta/Pz;
    syms R t
    radiuspart='1/(alphasigma/R-beta*Pz)';
    radiuspart=subs(radiuspart,'alphasigma',alphasigma);
    radiuspart=subs(radiuspart,'Pz',Pz);
    radiuspart=subs(radiuspart,'beta',beta');
    Ri=[linspace(R0,Rlim-0.00001,50) linspace(Rlim-0.00001,Rlim-eps,5)];
    for ni=1:length(Ri)
        f=int(radiuspart,'R',R0,Ri(ni));
        ti(ni)=eval(f)/M;
    end

    hold on
    plot(ti,Ri)
    axis([0 1.5e4 18 32])
end




%% ------- 3D CLOSE LOOK --------
%% Grain Vol
clear
figure; hold on
savedir='/home/cenna/Results/3DClose/curve_test_100/';
scale=1;
Mr=[50];
steps=[10];
tol=1*[1];
MPz=[0];
L=1;
m=1;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
alpha=2;
for n=1:size(Mr,2)
   
    voldata=importdata([savedir 'GrainVol' '.txt']);
    voldata(1:50,:)=[];
    Pz=voldata(:,3);
    movinginds=find((Pz==Pz(end))==1);
    time=voldata(movinginds,1);
    vol=voldata(movinginds,2);
    radius=(vol/pi/4*3).^(1/3);
    % radiusi=linspace(min(radius),max(radius),1000);
    % timei=interp1(time,radius,radiusi);
    % plot(time,radius)
    spline1 = spaps(time(1:steps(n):end),radius(1:steps(n):end),tol(n),3);
    Mvel=fnval(fnder(spline1,1),time(1:steps(n):end));
    % xlabel('Time')
    % ylabel('radius of circular grain')
    % figure
    Mvel=-Mvel;
    Mcurvature=alpha./radius(1:steps(n):end);
    plot(Mcurvature,Mvel/mobility/intenergy,'.')
    xlabel('\Delta G/ \sigma_{gb}');
    ylabel('V / (M \sigma_{gb})');
    grid on
    MPz=[MPz Pz(end)];
end
Normalized_Pz=MPz/intenergy
axis([0 0.2 0 0.2]);box on
hold on
plot([0 0.5],[0 0.5],'r')
title(savedir)

%% 
figure
plot(time(1:steps(n):end), Mcurvature,'.')
ylabel('\Delta G/ \sigma_{gb}')
xlabel('time')
%% plot quality of time-radius spline
figure
plot(time,radius)
xlabel('time')
ylabel('Radius')
hold on
spline1 = spaps(time(1:steps(n):end),radius(1:steps(n):end),tol(n),3);
plot(time(1:steps(n):end),fnval(spline1,time(1:steps(n):end)),'r-')
grid on
%% from r^2-r0^2 = kt
figure
plot(time,radius.^2)
xlabel('time')
ylabel('Radius ^2')

pp=polyfit(time, radius.^2,1)
slope_from_formula=-2*alpha*mobility*intenergy

%% analytical exprression for solute drag of flat interface
L=1;
m=2;
kappa=4;
a=0.25;
b=30;

mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);

vi=linspace(0,2,1000);
Ps=a*vi./(1+b*vi.^2);

hold on
plot(vi,Ps)
xlabel('V');ylabel('P_s')

%% ------------------------- PARTICLES N 3D -------------------------
%% particlesN 3D - Vollog
clear
ni=0;
L=1;
m=2;
kappa=4;
mobility=3/2*L*sqrt(2*kappa/m);
%mobility=2.77
intenergy=1/3*sqrt(2*m*kappa);
% fv=0.0467;
delx=0.5;
mboxsize=300;
% figure
DelGi=[ 500 400 300 200 100 60 50]
for ni=1:length(DelGi)
savedir=['/home/magnetadmin/Documents/Results/Np_Solute/Np300_m2_k4_r5_fv50_DelG' num2str(DelGi(ni)) '/'];
DelG=DelGi(ni)/1000;
nboxsize=mboxsize;
voldata=importdata([savedir 'Vollog' '.log']);
time=voldata(200:end,1);
vol=voldata(200:end,2);
fv=voldata(1,1)/1
vol=vol/(1-fv);
intlength=(vol-vol(1))/(mboxsize*nboxsize*delx*delx);

% hold on
% plot(time,intlength)
% xlabel('Time')
% ylabel('Position of the interface')

% finding Pz
pp=polyfit(time,intlength,1);
Mvel(ni)=-pp(1);
Pz(ni)=(mobility*DelG-Mvel(ni))/mobility;
end
figure
plot(DelGi/1000,Mvel,'o')
xlabel('\Delta G'); ylabel('Interface Velocity')
hold on; plot([0 max(DelGi)/1000],[0 mobility*max(DelGi)/1000],'r')

figure
plot(Mvel,Pz,'s')
xlabel('v','FontSize',14);ylabel('P_z','FontSize',14)

%% particlesN 3D - Energy Log
savedir='/home/magnetadmin/Documents/Results/Np_Solute/P1_m2_k4_r5/200/';

engdata=importdata([savedir 'Englog' '.log']);
time=engdata(:,1);
E=engdata(:,2);
figure; 
hold on
plot(time,E)
xlabel('Time')
ylabel('Total Energy of system')

%% particlesN 3D - Vol Log
clear
savedir='/home/magnetadmin/Documents/Results/Np_Solute/Np300_m2_k4_r5_fv0.05/200/';
delx=0.5;
mboxsize=300;
voldata=importdata([savedir 'Vollog' '.log']);
time=voldata(2:end,1);
vol=voldata(2:end,2);
% fv=voldata(1,1)/2
fv=0.0
vol=vol/(1-fv);
intlength=(vol)/(mboxsize*mboxsize*delx*delx);

  figure; 
hold on
plot(time,intlength,'r')
xlabel('Time')
ylabel('Interface Position')
%% finding Pz
DelG=0.20
pp=polyfit(time,intlength,1);
Mvel=-pp(1)
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
VwithNoPinning=mobility*DelG
Pz=(mobility*DelG-Mvel)/mobility


%% Serial Analysis for Interface_Particle Well
clear
ni=0;
DelG=[100]/1000;
for var=DelG
    ni=ni+1;
    savedir=['/home/cenna/Results/test/P50_m1_k2_r5' '/'];
    L=1;
    m=1;
    kappa=2;  
    mobility=3/2*L*sqrt(2*kappa/m);
    intenergy=1/3*sqrt(2*m*kappa);
    delx=0.5;
    mboxsize=70;
    nboxsize=mboxsize;
    voldata=importdata([savedir 'Vollog' '.log']);
    time=voldata(2:end,1);
    vol=voldata(2:end,2);
    intlength=vol/(mboxsize*nboxsize*delx*delx);
    engdata=importdata([savedir 'Englog' '.log']);
    time=engdata(:,1);
    E=engdata(:,2);
    intlength=mboxsize*delx-intlength;
    plot(intlength,E)
    hold on; grid on;
    xlabel('Interface Position')
    ylabel('Energy of the System')
    
    [Emax,imax]=extrema(E);

end

%% measuring error for no particle
figure
subplot(2,1,1)
plot(DelG,Mvel,'o')
hold on
plot([0 .1], mobility*[0 .1],'r')
xlabel('Driving Pressure')
ylabel('Interface Velocity')
subplot(2,1,2)
err=(mobility*DelG-Mvel)./(mobility*DelG)
stem(DelG,err)
ylabel('Relative Error')
