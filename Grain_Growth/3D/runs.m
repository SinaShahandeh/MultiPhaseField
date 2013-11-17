%% Topology Distributions
clear
Msavedir={['/media/Disk2/Results/3D/Fric300_Pz0_m1_k2_init0/']};
Mmboxsize=[300];
Mstart=[2000] ;
Msteps=[2000];
Mends=[54000];
for simi=1
    mboxsize=Mmboxsize(simi);
    nboxsize=mboxsize;
    lboxsize=mboxsize;

    delx=2;
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