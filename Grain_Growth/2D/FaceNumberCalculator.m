%% Face numbers analysis
function FaceNumberCalculator(savedir,start,steps,ends)
% clear
% savedir='/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz0/';
mkdir([savedir '/facenums'])
% steps=2000;
% start=2000;
% ends=428000;
Mtn=[start:steps:ends];
%  Mtn=[1000];
se=strel('square',3);
for tni=1:size(Mtn,2)
    tn=Mtn(tni)
    Inds=importdata([savedir 'Inds_' num2str(tn) '.txt']);
    mboxsize=size(Inds,1);nboxsize=size(Inds,2);
    % finding index of grains with non zero size:
    grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    [nonzerogi]=find(grains>10);
    % obtaining statistics
    facenum=zeros(size(nonzerogi));
    for gi=1:length(nonzerogi)
        ni=nonzerogi(gi);
        % search for ni in the Inds matrix
        BWInds=logical(zeros(size(Inds)));
        BWInds(Inds==ni)=1;
        [row,col]=find(BWInds==true);
        row(row<=2)=3; row(row>=mboxsize-2)=mboxsize-2;
        col(col<=2)=3; col(col>=nboxsize-2)=nboxsize-2;
        BWInds=BWInds(min(row)-2:max(row)+2,min(col)-2:max(col)+2);
        Inds2BW=Inds(min(row)-2:max(row)+2,min(col)-2:max(col)+2);
        BWInds=imdilate(BWInds,se);
        Inds2BW=Inds2BW.*double(BWInds);
        listi=Inds2BW(Inds2BW~=0);
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