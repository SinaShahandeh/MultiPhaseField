 %% Raw curvature histograms
clear
savedir='/home/magnetadmin/Documents/Results/3D/Fric300_Pz0.010_m1_k2_init4000/';
tn=900

Pz=0.013
L=1;
m=1;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m)
intenergy=1/3*sqrt(2*m*kappa)
tni=0;
start=5000;
step=5000;
ending=150000;
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
    
%     hist(log(Khist),50)
%     title([ 'Time step = ' num2str(tn) ], 'FontSize',16);h=get(gcf,'Children'); set(h(1),'LineWidth',1.5);set(h,'FontSize',14)
%     xlabel('Ln(Curvarure), Ln(K)','FontSize',16); ylabel('Length Fraction','FontSize',16);% pause(1)
%     hold on
%     plot(log([meanK(tni) meanK(tni)]), [0 1000],'g-' ,'LineWidth',2)
%     plot(log([Pz/intenergy Pz/intenergy]), [0 1000],'r-' ,'LineWidth',2)
%     hold off
%     mkdir([savedir 'curvarure_dist/'])
%     print([savedir 'curvarure_dist/' num2str(tni) '.png'],'-dpng','-r200',gcf)
 end

 bigly bigly
 
 %%

clear
savedir='/home/magnetadmin/Documents/Results/3D/Fric300_Pz0_m1_k2_init0_run2/';
start=5000;
step=5000;
ending=40000;
mboxsize=300;
nboxsize=300;
kboxsize=300;
m=1;
kappa=2;
a=sqrt(m/2/kappa);
delx=1;
delt=0.1;
pmax=8;
tni=0;
for tn=[start:step:ending]
    tni=tni+1;
    tn
    Kdata=importdata([savedir 'K_' num2str(tn) '.txt']);
    Kmat=zeros(mboxsize,nboxsize,nboxsize)*nan;
    for i=1:length(Kdata)
        if Kdata(i,5)>0  % select the positive component
            Kmat(Kdata(i,1)+1,Kdata(i,2)+1,Kdata(i,3)+1)=Kdata(i,5);
        end
    end
    Kmat(Kmat>0.06)=nan; % to detach segments at the triple junctions with high curvature
    Khist=[];
    Lseg=bwlabeln(~isnan(Kmat));
    segnum(tni)=0;
    for i=1:max(max(max(Lseg)))
        Segi=find(Lseg==i);
        if length(Segi)>1 % removing very small segments
            Ksegmat=Kmat(Segi);
            segnum(tni)=segnum(tni)+1;
            [f,ki]=hist((Kmat(Segi)),200); % finding position of maximum in histogram distribution and repeating the data for the whole segment
            [a, fmaxi]=max(f);
            rangek=Ksegmat(Ksegmat<ki(fmaxi)+0.01 & Ksegmat>ki(fmaxi)-0.01);
            Ksegi=mean(rangek)+zeros(length(Segi),1);
            Khist=[Khist Ksegi'];
        end
    end
    meanK(tni)=mean(Khist);
%     figure
%     hist((Khist),100)
    save( [savedir 'Khist_' num2str(tn) '.mat'], 'Khist')
end

figure
plot([start:step:ending]*delt, 1./meanK,'-o')
