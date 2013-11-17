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