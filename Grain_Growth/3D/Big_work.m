
%% 3D slice full reading (E)
clear
tn=1000;

savedir='/home/magnetadmin/Documents/Results/3D/Big/90/';
phidata=importdata([savedir 'Energy_' num2str(tn) '.txt']);
mboxsize=size(phidata,2);nboxsize=mboxsize;lboxsize=mboxsize;
% x=[0:nboxsize-1]*delx;
% y=[0:mboxsize-1]*delx;
% z=[0:lboxsize-1]*delx;
% [x,y,z]=meshgrid(x,y,z);

E=zeros(mboxsize,nboxsize,lboxsize);
for m=1:mboxsize
    for n=1:nboxsize
        for l=1:lboxsize
            E(m,n,l)=phidata(m+(l-1)*mboxsize,n);
        end
    end
end
clear phidata;
E(E>0.5)=nan;
% figure; slice(E,[1 (mboxsize)],[1 (nboxsize)],[1 (lboxsize)]);colormap gray; shading flat; axis equal; box on; colorbar; title ([savedir ' (tn=' num2str(tn) ')'])
%  axis([1 (mboxsize) 1 (nboxsize) 1 (lboxsize)])
sliceomatic(E)


%% Grain Size from GrainStat files

clear
n=1;
for simnum=[10 20 30 40 50 60 70 80 90]
savedir=['/home/magnetadmin/Documents/Results/3D/Big/' num2str(simnum) '/'];
delt=0.1;
delx=1;
start=100; steps=100;
Mtn=[start:steps:1300000];
clear GrainStat Volbig Volbar grains 
try
    for tni=1:size(Mtn,2)
        tn=Mtn(tni);
        GrainStat{tni}=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
        grains=GrainStat{tni};
        grains(grains<3^3)=[];
        volbig(tni)=max(grains)*delx^3;
        grains(grains==max(grains))=[];
        volbar(tni)=mean(grains)*delx^3;
        numgrais(tni)=length(grains);
    end
catch
end     
Mtn=[start:steps:tn-steps];
timevec=Mtn*delt;
%figure
%plot(timevec,volbar,'o')
%ylabel('Mean Volume of Grains');xlabel('Time')
R_bar=((3/4/pi)^(1/3))*volbar.^(1/3);
R_big=((3/4/pi)^(1/3))*volbig.^(1/3);
% figure
% hold on
% plot(timevec,R_big./R_bar,'ro')
% ylabel('(Equivalent Radius)');xlabel('Time');grid on
pp=polyfit(timevec(8:end),R_big(8:end),1)
slope(n) = pp(1);
initialsize(n)= R_big(1);
finalsize(n)=R_big(end);
R_bar_end(n)=R_bar(end);
ratio(n)=R_big(end)/R_bar(end);
n=n+1;
end
%%
figure
plot(ratio,slope,'s','LineWidth',1)
grid on
ylabel('Growth Speed of Large grain, (dR_A/dt)','FontSize',14)
xlabel('R_A / R_m','FontSize',14)
set(get(gcf,'Children'),'LineWidth',1.5)
set(get(gcf,'Children'),'FontSize',12)


%% total energy
TotalE=[2.42394e+06 2.40656e+06 2.37874e+06 2.33828e+06 2.28525e+06 2.20938e+06 2.11101e+06 1.99848e+06 1.85503e+06]

figure
plot(finalsize,TotalE,'s','LineWidth',1)
grid on
ylabel('Total Energy of the System, (E)','FontSize',14)
xlabel('R_A','FontSize',14)
set(get(gcf,'Children'),'LineWidth',1.5)
set(get(gcf,'Children'),'FontSize',12)

L=1;
m=1;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa)

em=intenergy/2/mean(R_bar_end)


%% Finding neighbors of a grain (number of faces)
clear
 nii=1;
for simnum=[10]
    mboxsize=300;
    nboxsize=mboxsize;
    lboxsize=mboxsize;
    nuclein=abs(mboxsize*nboxsize*lboxsize/1000);
    savedir=['/home/magnetadmin/Documents/Results/3D/Big/' num2str(simnum) '/'];
    tn=1000;
    % finding the bigest grain index:
    grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    [a,ind]=max(grains);grains(ind)=0;
    se=strel(ones(2,2,2));
    Inds=importdata([savedir 'Inds_' num2str(tn) '.txt']);
    BWInds=uint8(zeros(size(Inds)));
    ni=ind
    % search for ni in the Inds matrix
    BWInds(Inds==ni)=1;
    BWInds=imdilate(BWInds,se);
    Inds2BW=Inds.*double(BWInds);
    listi=Inds(Inds2BW~=0);
    facenum(nii)=0
    for i=1:length(listi)
        indi=listi(i);
        if indi~=0
            facenum(nii)=facenum(nii)+1;
            listi(listi==indi)=0;
        end
    end
    nii=nii+1;
end


%%
facenum=[83   133   186   267   362   446   570   689    819]
figure
plot(ratio,facenum,'s','LineWidth',1)
grid on
xlabel('Ratio of R_A/R_m','FontSize',14)
ylabel('Number of faces of large grain','FontSize',14)
set(get(gcf,'Children'),'LineWidth',1.5)
set(get(gcf,'Children'),'FontSize',12)

figure
plot(facenum,finalsize,'s','LineWidth',1)
grid on
ylabel('Large Grain Equivalent Radius','FontSize',14)
xlabel('Number of faces of large grain','FontSize',14)
set(get(gcf,'Children'),'LineWidth',1.5)
set(get(gcf,'Children'),'FontSize',12)


%% Topology Distributions
clear
for simnum=[10 20 30 40 50 60 70 80 90]
    mboxsize=300;
    nboxsize=mboxsize;
    lboxsize=mboxsize;
    nuclein=abs(mboxsize*nboxsize*lboxsize/1000);
    delx=1;
    delt=0.1;
    savedir=['/home/magnetadmin/Documents/Results/3D/Big/' num2str(simnum) '/'];
    mkdir([savedir '/facenums'])
    se=strel(ones(2,2,2));
    tn=1000;
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


