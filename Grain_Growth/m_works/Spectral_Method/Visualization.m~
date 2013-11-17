%% Visualization

%% 2D cross sections
savedir='/home/cenna/CMPE/Results/grgr3d_1/'
for tn=1:1000
    load([savedir num2str(tn) '.mat'])
   pcolor(sumetasqu(:,:,1));
            caxis([0 1]);
            colormap gray
            shading('interp')
            axis off
            box off 
            set(gca,'DataAspectratio',[1 1 1])
            pause(0.01);
end

%% 3D slice
clear
savedir='/media/disk/sim_res/grgr3d_1/'
for tn=[10:10:1000]
    load([savedir num2str(tn) '.mat'])
    slice(sumetasqu,[1 size(sumetasqu,1)],[1 size(sumetasqu,2)],[1 size(sumetasqu,3)])
    colormap gray
    shading interp
    axis equal
    box on
    title (['timestep=' num2str(tn)])
    print(['/media/disk/sim_res/grgr3d_1/' num2str(tn) '.png']...
        ,'-dpng','-r200','-f1')
end
