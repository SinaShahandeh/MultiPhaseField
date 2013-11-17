function GrainR_t_3D(savedir,start,steps,ends,delt)

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
        GrainStat=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
        grains=GrainStat;
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
timevec=timevec-timevec(1);
figure
hold on
plot(timevec,Rbar,'b.')
ylabel('Equivalent Radius');xlabel('Time');grid on
% 
title(savedir);


set(gca,'Box', 'on')
set(gca,'FontSize',13,'LineWidth',2)

Rinit=Rbar(1)
