%% over time steps
clear
savedir='/home/magnetadmin/Documents/Results/';
start=5000;
step=5000;
ending=50000;
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
    tni=tni+1
    Kdata=importdata([savedir 'K_' num2str(tn) '.txt']);
    Kmat=zeros(mboxsize,nboxsize,nboxsize)*nan;
    for i=1:length(Kdata)
        if Kdata(i,5)>0
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
            [f,ki]=hist((Kmat(Segi)),500); % finding position of maximum in histogram distribution and repeating the data for the whole segment
            [a, fmaxi]=max(f);
            rangek=Ksegmat(Ksegmat<ki(fmaxi)+0.005 & Ksegmat>ki(fmaxi)-0.005);
            Ksegi=mean(rangek)+zeros(length(Segi),1);
            Khist=[Khist Ksegi'];
        end
    end
    meanK(tni)=mean(Khist);
    figure
    hist((Khist),50)
    save( [savedir 'Khist0.02_' num2str(tn) '.mat'], 'Khist')
end