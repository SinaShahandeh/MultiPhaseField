clear
figure
savedir='/media/disk/sim_res/DomeDrivingForce_solute_3/';
Mparam=linspace(50,500,10);
for filenum=1:length(Mparam)
    param=Mparam(filenum)
    [InterfaceVel]=realtime_particle_friction_solute(param,filenum,savedir);
    MInterfaceVel(filenum)=InterfaceVel;
%     Mpp(filenum)=pp
end
figure
plot(Mparam,MInterfaceVel)
grid

