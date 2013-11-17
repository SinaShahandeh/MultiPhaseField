%% === Friction of fine particles on the interface motion ===
%% changing Pz and track speed of the interface
clear
figure
savedir='/media/disk/sim_res/DomeDrivingForce_particle_03/';
mkdir(savedir)
% savedir= [pwd '/']
Mparam=linspace(10,50,20);
for filenum=1:length(Mparam)
    param=Mparam(filenum)
    [InterfaceVel]=realtime_particle_friction(param,filenum,savedir,0.05);
    MInterfaceVel(filenum)=InterfaceVel;
%     Mpp(filenum)=pp
end
figure
plot(Mparam,MInterfaceVel)
grid
save([savedir 'full.mat'])