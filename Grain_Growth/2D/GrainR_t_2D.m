function GrainR_t_2D(savedir,start,steps,ends,delt)

figure

delx=1;

L=1;
m=2;
kappa=4;
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

% figure
plot(timevec,Rbar,'.')
ylabel('Mean Radius of Grains');xlabel('Time')
%  figure
 hold on
% plot(tau,Rbar.^2,'bd','LineWidth',1.5)
% ylabel('<R>^2','FontSize',14);xlabel('\tau','FontSize',14);grid on;
% hold on

set(gca,'Box', 'on')
set(gca,'FontSize',13,'LineWidth',2)

Rinit=Rbar(1)