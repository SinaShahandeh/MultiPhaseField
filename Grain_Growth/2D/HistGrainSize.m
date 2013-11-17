%% size histograms from grain stats
function Rbar=HistGrainSize(savedir,tn)
% clear
% savedir='/home/cenna/Results/2D/Fric1000s2_m2_k3_phi_init1000/20_2/';
delt=0.1
delx=1;
% start=254000
% steps=2000
% Mtn=[start:steps:2000000];
figure
% tni=1;tn=10000
grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
grains(grains<10)=[];
numgrains=length(grains)
R=(grains/pi).^0.5;
Rbar=mean(R)*delx;
% logarithmic scale
Ri=linspace(log(5)/log(Rbar),log((max(R)+5))/log(Rbar),30);
N=hist(log(R)/log(Rbar),Ri);
bar(Ri,N)



title(['Time= ' num2str(tn)])

