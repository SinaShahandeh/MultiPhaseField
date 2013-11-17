clear
savedir='/media/Disk2/Results/2D/Fric2000_m2_k3_init2000/30/';
start=10000;
step=10000;
ending=600000;
mboxsize=2000;
nboxsize=2000;
m=2;
kappa=3;
a=sqrt(m/2/kappa);
delx=1;
delt=0.1;
pmax=5;
eta=zeros(mboxsize,nboxsize,pmax);
inds=eta;
tni=0;
for tn=[start:step:ending]
    tni=tni+1
    for p=1:pmax
        eta(:,:,p)=importdata([savedir 'Eta_' num2str(p-1) '_' num2str(tn) '.txt']);
        inds(:,:,p)=importdata([savedir 'inds_' num2str(p-1) '_' num2str(tn) '.txt']);
    end
    GrainStat=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    BiggestGrainD=2*sqrt(max(GrainStat)/pi); % to repeat the the system at the boundaries to include the largest grain if it falls on the boundary
    BigGrains=find(GrainStat>100)';
    gni=0; % counter for grains
    for gni=1:length(BigGrains)
        gi=BigGrains(gni);
        mboxsize=size(eta,1);nboxsize=size(eta,2);
        etani=zeros(mboxsize,nboxsize);
        giindex=find(inds==gi); % find grain gi in the inds, giindex is all the pixels in the inds that has gi
        for ii=1:length(giindex) % create etani matrix to contain only grain gi
            p=fix(giindex(ii)/mboxsize/nboxsize-eps)+1;
            j=fix((giindex(ii)-(p-1)*mboxsize*nboxsize)/nboxsize-eps)+1;
            i=giindex(ii)-(p-1)*mboxsize*nboxsize-(j-1)*nboxsize+1;
            etani(i,j)=eta(giindex(ii));
        end
        [in,jn]=find(etani>0.1);
        jmin=min(jn);jmax=max(jn);imin=min(in);imax=max(in);
        % Transfering the domain to make the bounudary grains at the centre
%         if jmax+BiggestGrainD>nboxsize
%             etani=[etani(:,nboxsize/2:end)  etani(:,1:nboxsize/2-1)];
%         end
%         if jmin-BiggestGrainD<1
%             etani=[etani(:,nboxsize/2:end)  etani(:,1:nboxsize/2-1)];
%         end
%         if imax+BiggestGrainD>mboxsize
%             etani=[etani(mboxsize/2:end,:); etani(1:mboxsize/2-1,:)];
%         end
%         if imin-BiggestGrainD<1
%             etani=[etani(mboxsize/2:end,:); etani(1:mboxsize/2-1,:)];
%         end
        if (imax==mboxsize) && (imin==1)
            etani=[etani(mboxsize/2:end,:); etani(1:mboxsize/2-1,:)];
        end
        if (jmin==1) && (jmax==nboxsize)
            etani=[etani(:,nboxsize/2:end)  etani(:,1:nboxsize/2-1)];
        end
        % cropping the grain from the whole domain
        [in,jn]=find(etani>0.1);
        jmin=min(jn);jmax=max(jn);imin=min(in);imax=max(in);
        etani=etani(imin:imax,jmin:jmax);

        mboxsize=size(etani,1); nboxsize=size(etani,2);
        %% calculating curvaturte from etani
        del2eta_1=zeros(mboxsize,nboxsize);
        for i=1:mboxsize
            for j=1:nboxsize
                %              del2eta(i,j)=1/delx^2*(etani(indg(i+1,mboxsize),j)-2*etani(i,j)+etani(indg(i-1,mboxsize),j))...
                %                +1/delx^2*(etani(i,indg(j+1,nboxsize))-2*etani(i,j)+etani(i,indg(j-1,nboxsize)));
                sumeta=etani(indg(i+1,mboxsize),j)+etani(indg(i-1,mboxsize),j)+etani(i,indg(j+1,nboxsize))+etani(i,indg(j-1,nboxsize));
                sumeta2=etani(indg(i+1,mboxsize),indg(j+1,nboxsize))+etani(indg(i-1,mboxsize),indg(j+1,nboxsize))+etani(indg(i+1,mboxsize),indg(j-1,nboxsize))+etani(indg(i-1,mboxsize),indg(j-1,nboxsize));
                del2eta_1(i,j)=2/3/(delx*delx)*(sumeta+0.25*sumeta2-5*etani(i,j));
            end
        end
        del2eta=1*del2eta_1;
        gradeta=-0.5*sqrt(m/2/kappa)*(1-(2*etani-1).^2);
        grad2eta=-sqrt(m/2/kappa)^2*(2*etani-1).*(1-(2*etani-1).^2);
        H=(del2eta-grad2eta)./gradeta;
        H(abs(H)>0.07)=nan;H(etani>0.8)=nan;H(etani<0.2)=nan; % filtering H
        %% ploting
%         gcf;clf;
%         surf(H)
%         view([0 90])
%         shading flat
%         colorbar;axis equal; box on
%         xlabel('X');ylabel('Y'); title([  'Grain number ' num2str(gi)])
%         ax=get(gcf,'children');set(ax(2),'Clim',[-0.03 0.03])
%         pause(0.1)
%         MH(:,:,gni)=H;
        %% sampling based on middle poi2nt of the segment
        BWsegs=~isnan(H);
        Lsegs=bwlabel(BWsegs,8); % 8 neighbor
        prop=regionprops(Lsegs,'Area','BoundingBox');
        Hsamples=[];
        MKSegs=zeros(1,length(prop));
        for segi=1:length(prop)
            % find minimum of H in the segment
            Hseg=H(fix(prop(segi).BoundingBox(2))+1:fix(prop(segi).BoundingBox(2))+fix(prop(segi).BoundingBox(4)), fix(prop(segi).BoundingBox(1))+1:fix(prop(segi).BoundingBox(1))+fix(prop(segi).BoundingBox(3)));
            [Hsegmin]=min(min(Hseg));
            [Hsegmini,Hsegminj]=find(Hseg==Hsegmin);
            selbnd=[fix(0.1*prop(segi).BoundingBox(4)) fix(0.1*prop(segi).BoundingBox(3))]+1; % data selection bounding box is 0.2 of segment bounding box
            try
            Ksegi=Hseg(Hsegmini-selbnd(1):Hsegmini+selbnd(1), Hsegminj-selbnd(2):Hsegminj+selbnd(2));
             nanelems=sum(sum(isnan(Ksegi)));
            Ksegi(isnan(Ksegi))=0;
            Ksegi=sum(sum(Ksegi))/(size(Ksegi,1)*size(Ksegi,2)-nanelems); % average of the Ksegi matrix without considering nan elements
            catch
                Ksegi=Hsegmin; % if the segment is too small then take the min for the whole segment
            end
            if isnan(Ksegi)
                disp('Warning: A too-much-urved interface!!!')
            end
            Hsamples=[Hsamples ones(1,prop(segi).Area)*Ksegi]; % repeats the Ksegi as much as the area of the segment
            MKSegs(segi)=Ksegi;
        end
        
         %%  sampling of H from the whole structure
%         [i,j]=find(~isnan(H)); % converting the H matrix with lots of NaNs to a vercor (Hsamples) for just one grain
%         Hsamples=zeros(1,length(i));
%         for n=1:length(i)
%             Hsamples(n)=H(i(n),j(n));
%         end
        MHsamples{gni}=Hsamples; % keeps all the grains data
    end
    
    %% analysing MHSamples (curvature histograms)
    AllHsamples=[]; % for all the structure conbining all inviditual grains curvarue data to Hsamples
    for ii=1:length(MHsamples)
        AllHsamples=[AllHsamples MHsamples{ii}];
    end
    %     figure
    %     hist (AllHsamples,100)
    Mtn_K{tni}=AllHsamples;
    MHsamples={};
end
clear eta inds
save([ savedir 'curvature_data.mat'])

%%


clear
savedir='/media/Disk2/Results/2D/Fric2000_m2_k3_init2000/20/';
start=10000;
step=10000;
ending=600000;
mboxsize=2000;
nboxsize=2000;
m=2;
kappa=3;
a=sqrt(m/2/kappa);
delx=1;
delt=0.1;
pmax=5;
eta=zeros(mboxsize,nboxsize,pmax);
inds=eta;
tni=0;
for tn=[start:step:ending]
    tni=tni+1
    for p=1:pmax
        eta(:,:,p)=importdata([savedir 'Eta_' num2str(p-1) '_' num2str(tn) '.txt']);
        inds(:,:,p)=importdata([savedir 'inds_' num2str(p-1) '_' num2str(tn) '.txt']);
    end
    GrainStat=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    BiggestGrainD=2*sqrt(max(GrainStat)/pi); % to repeat the the system at the boundaries to include the largest grain if it falls on the boundary
    BigGrains=find(GrainStat>100)';
    gni=0; % counter for grains
    for gni=1:length(BigGrains)
        gi=BigGrains(gni);
        mboxsize=size(eta,1);nboxsize=size(eta,2);
        etani=zeros(mboxsize,nboxsize);
        giindex=find(inds==gi); % find grain gi in the inds, giindex is all the pixels in the inds that has gi
        for ii=1:length(giindex) % create etani matrix to contain only grain gi
            p=fix(giindex(ii)/mboxsize/nboxsize-eps)+1;
            j=fix((giindex(ii)-(p-1)*mboxsize*nboxsize)/nboxsize-eps)+1;
            i=giindex(ii)-(p-1)*mboxsize*nboxsize-(j-1)*nboxsize+1;
            etani(i,j)=eta(giindex(ii));
        end
        [in,jn]=find(etani>0.1);
        jmin=min(jn);jmax=max(jn);imin=min(in);imax=max(in);
        % Transfering the domain to make the bounudary grains at the centre
%         if jmax+BiggestGrainD>nboxsize
%             etani=[etani(:,nboxsize/2:end)  etani(:,1:nboxsize/2-1)];
%         end
%         if jmin-BiggestGrainD<1
%             etani=[etani(:,nboxsize/2:end)  etani(:,1:nboxsize/2-1)];
%         end
%         if imax+BiggestGrainD>mboxsize
%             etani=[etani(mboxsize/2:end,:); etani(1:mboxsize/2-1,:)];
%         end
%         if imin-BiggestGrainD<1
%             etani=[etani(mboxsize/2:end,:); etani(1:mboxsize/2-1,:)];
%         end
        if (imax==mboxsize) && (imin==1)
            etani=[etani(mboxsize/2:end,:); etani(1:mboxsize/2-1,:)];
        end
        if (jmin==1) && (jmax==nboxsize)
            etani=[etani(:,nboxsize/2:end)  etani(:,1:nboxsize/2-1)];
        end
        % cropping the grain from the whole domain
        [in,jn]=find(etani>0.1);
        jmin=min(jn);jmax=max(jn);imin=min(in);imax=max(in);
        etani=etani(imin:imax,jmin:jmax);

        mboxsize=size(etani,1); nboxsize=size(etani,2);
        %% calculating curvaturte from etani
        del2eta_1=zeros(mboxsize,nboxsize);
        for i=1:mboxsize
            for j=1:nboxsize
                %              del2eta(i,j)=1/delx^2*(etani(indg(i+1,mboxsize),j)-2*etani(i,j)+etani(indg(i-1,mboxsize),j))...
                %                +1/delx^2*(etani(i,indg(j+1,nboxsize))-2*etani(i,j)+etani(i,indg(j-1,nboxsize)));
                sumeta=etani(indg(i+1,mboxsize),j)+etani(indg(i-1,mboxsize),j)+etani(i,indg(j+1,nboxsize))+etani(i,indg(j-1,nboxsize));
                sumeta2=etani(indg(i+1,mboxsize),indg(j+1,nboxsize))+etani(indg(i-1,mboxsize),indg(j+1,nboxsize))+etani(indg(i+1,mboxsize),indg(j-1,nboxsize))+etani(indg(i-1,mboxsize),indg(j-1,nboxsize));
                del2eta_1(i,j)=2/3/(delx*delx)*(sumeta+0.25*sumeta2-5*etani(i,j));
            end
        end
        del2eta=1*del2eta_1;
        gradeta=-0.5*sqrt(m/2/kappa)*(1-(2*etani-1).^2);
        grad2eta=-sqrt(m/2/kappa)^2*(2*etani-1).*(1-(2*etani-1).^2);
        H=(del2eta-grad2eta)./gradeta;
        H(abs(H)>0.07)=nan;H(etani>0.8)=nan;H(etani<0.2)=nan; % filtering H
        %% ploting
%         gcf;clf;
%         surf(H)
%         view([0 90])
%         shading flat
%         colorbar;axis equal; box on
%         xlabel('X');ylabel('Y'); title([  'Grain number ' num2str(gi)])
%         ax=get(gcf,'children');set(ax(2),'Clim',[-0.03 0.03])
%         pause(0.1)
%         MH(:,:,gni)=H;
        %% sampling based on middle poi2nt of the segment
        BWsegs=~isnan(H);
        Lsegs=bwlabel(BWsegs,8); % 8 neighbor
        prop=regionprops(Lsegs,'Area','BoundingBox');
        Hsamples=[];
        MKSegs=zeros(1,length(prop));
        for segi=1:length(prop)
            % find minimum of H in the segment
            Hseg=H(fix(prop(segi).BoundingBox(2))+1:fix(prop(segi).BoundingBox(2))+fix(prop(segi).BoundingBox(4)), fix(prop(segi).BoundingBox(1))+1:fix(prop(segi).BoundingBox(1))+fix(prop(segi).BoundingBox(3)));
            [Hsegmin]=min(min(Hseg));
            [Hsegmini,Hsegminj]=find(Hseg==Hsegmin);
            selbnd=[fix(0.1*prop(segi).BoundingBox(4)) fix(0.1*prop(segi).BoundingBox(3))]+1; % data selection bounding box is 0.2 of segment bounding box
            try
            Ksegi=Hseg(Hsegmini-selbnd(1):Hsegmini+selbnd(1), Hsegminj-selbnd(2):Hsegminj+selbnd(2));
             nanelems=sum(sum(isnan(Ksegi)));
            Ksegi(isnan(Ksegi))=0;
            Ksegi=sum(sum(Ksegi))/(size(Ksegi,1)*size(Ksegi,2)-nanelems); % average of the Ksegi matrix without considering nan elements
            catch
                Ksegi=Hsegmin; % if the segment is too small then take the min for the whole segment
            end
            if isnan(Ksegi)
                disp('Warning: A too-much-urved interface!!!')
            end
            Hsamples=[Hsamples ones(1,prop(segi).Area)*Ksegi]; % repeats the Ksegi as much as the area of the segment
            MKSegs(segi)=Ksegi;
        end
        
         %%  sampling of H from the whole structure
%         [i,j]=find(~isnan(H)); % converting the H matrix with lots of NaNs to a vercor (Hsamples) for just one grain
%         Hsamples=zeros(1,length(i));
%         for n=1:length(i)
%             Hsamples(n)=H(i(n),j(n));
%         end
        MHsamples{gni}=Hsamples; % keeps all the grains data
    end
    
    %% analysing MHSamples (curvature histograms)
    AllHsamples=[]; % for all the structure conbining all inviditual grains curvarue data to Hsamples
    for ii=1:length(MHsamples)
        AllHsamples=[AllHsamples MHsamples{ii}];
    end
    %     figure
    %     hist (AllHsamples,100)
    Mtn_K{tni}=AllHsamples;
    MHsamples={};
end
clear eta inds
save([ savedir 'curvature_data.mat'])

%%



clear
savedir='/media/Disk2/Results/2D/Fric2000_m2_k3_init2000/15/';
start=10000;
step=10000;
ending=600000;
mboxsize=2000;
nboxsize=2000;
m=2;
kappa=3;
a=sqrt(m/2/kappa);
delx=1;
delt=0.1;
pmax=5;
eta=zeros(mboxsize,nboxsize,pmax);
inds=eta;
tni=0;
for tn=[start:step:ending]
    tni=tni+1
    for p=1:pmax
        eta(:,:,p)=importdata([savedir 'Eta_' num2str(p-1) '_' num2str(tn) '.txt']);
        inds(:,:,p)=importdata([savedir 'inds_' num2str(p-1) '_' num2str(tn) '.txt']);
    end
    GrainStat=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    BiggestGrainD=2*sqrt(max(GrainStat)/pi); % to repeat the the system at the boundaries to include the largest grain if it falls on the boundary
    BigGrains=find(GrainStat>100)';
    gni=0; % counter for grains
    for gni=1:length(BigGrains)
        gi=BigGrains(gni);
        mboxsize=size(eta,1);nboxsize=size(eta,2);
        etani=zeros(mboxsize,nboxsize);
        giindex=find(inds==gi); % find grain gi in the inds, giindex is all the pixels in the inds that has gi
        for ii=1:length(giindex) % create etani matrix to contain only grain gi
            p=fix(giindex(ii)/mboxsize/nboxsize-eps)+1;
            j=fix((giindex(ii)-(p-1)*mboxsize*nboxsize)/nboxsize-eps)+1;
            i=giindex(ii)-(p-1)*mboxsize*nboxsize-(j-1)*nboxsize+1;
            etani(i,j)=eta(giindex(ii));
        end
        [in,jn]=find(etani>0.1);
        jmin=min(jn);jmax=max(jn);imin=min(in);imax=max(in);
        % Transfering the domain to make the bounudary grains at the centre
%         if jmax+BiggestGrainD>nboxsize
%             etani=[etani(:,nboxsize/2:end)  etani(:,1:nboxsize/2-1)];
%         end
%         if jmin-BiggestGrainD<1
%             etani=[etani(:,nboxsize/2:end)  etani(:,1:nboxsize/2-1)];
%         end
%         if imax+BiggestGrainD>mboxsize
%             etani=[etani(mboxsize/2:end,:); etani(1:mboxsize/2-1,:)];
%         end
%         if imin-BiggestGrainD<1
%             etani=[etani(mboxsize/2:end,:); etani(1:mboxsize/2-1,:)];
%         end
        if (imax==mboxsize) && (imin==1)
            etani=[etani(mboxsize/2:end,:); etani(1:mboxsize/2-1,:)];
        end
        if (jmin==1) && (jmax==nboxsize)
            etani=[etani(:,nboxsize/2:end)  etani(:,1:nboxsize/2-1)];
        end
        % cropping the grain from the whole domain
        [in,jn]=find(etani>0.1);
        jmin=min(jn);jmax=max(jn);imin=min(in);imax=max(in);
        etani=etani(imin:imax,jmin:jmax);

        mboxsize=size(etani,1); nboxsize=size(etani,2);
        %% calculating curvaturte from etani
        del2eta_1=zeros(mboxsize,nboxsize);
        for i=1:mboxsize
            for j=1:nboxsize
                %              del2eta(i,j)=1/delx^2*(etani(indg(i+1,mboxsize),j)-2*etani(i,j)+etani(indg(i-1,mboxsize),j))...
                %                +1/delx^2*(etani(i,indg(j+1,nboxsize))-2*etani(i,j)+etani(i,indg(j-1,nboxsize)));
                sumeta=etani(indg(i+1,mboxsize),j)+etani(indg(i-1,mboxsize),j)+etani(i,indg(j+1,nboxsize))+etani(i,indg(j-1,nboxsize));
                sumeta2=etani(indg(i+1,mboxsize),indg(j+1,nboxsize))+etani(indg(i-1,mboxsize),indg(j+1,nboxsize))+etani(indg(i+1,mboxsize),indg(j-1,nboxsize))+etani(indg(i-1,mboxsize),indg(j-1,nboxsize));
                del2eta_1(i,j)=2/3/(delx*delx)*(sumeta+0.25*sumeta2-5*etani(i,j));
            end
        end
        del2eta=1*del2eta_1;
        gradeta=-0.5*sqrt(m/2/kappa)*(1-(2*etani-1).^2);
        grad2eta=-sqrt(m/2/kappa)^2*(2*etani-1).*(1-(2*etani-1).^2);
        H=(del2eta-grad2eta)./gradeta;
        H(abs(H)>0.07)=nan;H(etani>0.8)=nan;H(etani<0.2)=nan; % filtering H
        %% ploting
%         gcf;clf;
%         surf(H)
%         view([0 90])
%         shading flat
%         colorbar;axis equal; box on
%         xlabel('X');ylabel('Y'); title([  'Grain number ' num2str(gi)])
%         ax=get(gcf,'children');set(ax(2),'Clim',[-0.03 0.03])
%         pause(0.1)
%         MH(:,:,gni)=H;
        %% sampling based on middle poi2nt of the segment
        BWsegs=~isnan(H);
        Lsegs=bwlabel(BWsegs,8); % 8 neighbor
        prop=regionprops(Lsegs,'Area','BoundingBox');
        Hsamples=[];
        MKSegs=zeros(1,length(prop));
        for segi=1:length(prop)
            % find minimum of H in the segment
            Hseg=H(fix(prop(segi).BoundingBox(2))+1:fix(prop(segi).BoundingBox(2))+fix(prop(segi).BoundingBox(4)), fix(prop(segi).BoundingBox(1))+1:fix(prop(segi).BoundingBox(1))+fix(prop(segi).BoundingBox(3)));
            [Hsegmin]=min(min(Hseg));
            [Hsegmini,Hsegminj]=find(Hseg==Hsegmin);
            selbnd=[fix(0.1*prop(segi).BoundingBox(4)) fix(0.1*prop(segi).BoundingBox(3))]+1; % data selection bounding box is 0.2 of segment bounding box
            try
            Ksegi=Hseg(Hsegmini-selbnd(1):Hsegmini+selbnd(1), Hsegminj-selbnd(2):Hsegminj+selbnd(2));
             nanelems=sum(sum(isnan(Ksegi)));
            Ksegi(isnan(Ksegi))=0;
            Ksegi=sum(sum(Ksegi))/(size(Ksegi,1)*size(Ksegi,2)-nanelems); % average of the Ksegi matrix without considering nan elements
            catch
                Ksegi=Hsegmin; % if the segment is too small then take the min for the whole segment
            end
            if isnan(Ksegi)
                disp('Warning: A too-much-urved interface!!!')
            end
            Hsamples=[Hsamples ones(1,prop(segi).Area)*Ksegi]; % repeats the Ksegi as much as the area of the segment
            MKSegs(segi)=Ksegi;
        end
        
         %%  sampling of H from the whole structure
%         [i,j]=find(~isnan(H)); % converting the H matrix with lots of NaNs to a vercor (Hsamples) for just one grain
%         Hsamples=zeros(1,length(i));
%         for n=1:length(i)
%             Hsamples(n)=H(i(n),j(n));
%         end
        MHsamples{gni}=Hsamples; % keeps all the grains data
    end
    
    %% analysing MHSamples (curvature histograms)
    AllHsamples=[]; % for all the structure conbining all inviditual grains curvarue data to Hsamples
    for ii=1:length(MHsamples)
        AllHsamples=[AllHsamples MHsamples{ii}];
    end
    %     figure
    %     hist (AllHsamples,100)
    Mtn_K{tni}=AllHsamples;
    MHsamples={};
end
clear eta inds
save([ savedir 'curvature_data.mat'])

%%


clear
savedir='/home/cenna/Results/2D/Fric1000s2_m2_k3_phi_init1000/50/';
start=500;
step=500;
ending=32000;
mboxsize=2000;
nboxsize=2000;
m=2;
kappa=3;
a=sqrt(m/2/kappa);
delx=1;
delt=0.1;
pmax=5;
eta=zeros(mboxsize,nboxsize,pmax);
inds=eta;
tni=0;
for tn=[start:step:ending]
    tni=tni+1
    for p=1:pmax
        eta(:,:,p)=importdata([savedir 'Eta_' num2str(p-1) '_' num2str(tn) '.txt']);
        inds(:,:,p)=importdata([savedir 'inds_' num2str(p-1) '_' num2str(tn) '.txt']);
    end
    GrainStat=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
    BiggestGrainD=2*sqrt(max(GrainStat)/pi); % to repeat the the system at the boundaries to include the largest grain if it falls on the boundary
    BigGrains=find(GrainStat>100)';
    gni=0; % counter for grains
    for gni=1:length(BigGrains)
        gi=BigGrains(gni);
        mboxsize=size(eta,1);nboxsize=size(eta,2);
        etani=zeros(mboxsize,nboxsize);
        giindex=find(inds==gi); % find grain gi in the inds, giindex is all the pixels in the inds that has gi
        for ii=1:length(giindex) % create etani matrix to contain only grain gi
            p=fix(giindex(ii)/mboxsize/nboxsize-eps)+1;
            j=fix((giindex(ii)-(p-1)*mboxsize*nboxsize)/nboxsize-eps)+1;
            i=giindex(ii)-(p-1)*mboxsize*nboxsize-(j-1)*nboxsize+1;
            etani(i,j)=eta(giindex(ii));
        end
        [in,jn]=find(etani>0.1);
        jmin=min(jn);jmax=max(jn);imin=min(in);imax=max(in);
        % Transfering the domain to make the bounudary grains at the centre
%         if jmax+BiggestGrainD>nboxsize
%             etani=[etani(:,nboxsize/2:end)  etani(:,1:nboxsize/2-1)];
%         end
%         if jmin-BiggestGrainD<1
%             etani=[etani(:,nboxsize/2:end)  etani(:,1:nboxsize/2-1)];
%         end
%         if imax+BiggestGrainD>mboxsize
%             etani=[etani(mboxsize/2:end,:); etani(1:mboxsize/2-1,:)];
%         end
%         if imin-BiggestGrainD<1
%             etani=[etani(mboxsize/2:end,:); etani(1:mboxsize/2-1,:)];
%         end
        if (imax==mboxsize) && (imin==1)
            etani=[etani(mboxsize/2:end,:); etani(1:mboxsize/2-1,:)];
        end
        if (jmin==1) && (jmax==nboxsize)
            etani=[etani(:,nboxsize/2:end)  etani(:,1:nboxsize/2-1)];
        end
        % cropping the grain from the whole domain
        [in,jn]=find(etani>0.1);
        jmin=min(jn);jmax=max(jn);imin=min(in);imax=max(in);
        etani=etani(imin:imax,jmin:jmax);

        mboxsize=size(etani,1); nboxsize=size(etani,2);
        %% calculating curvaturte from etani
        del2eta_1=zeros(mboxsize,nboxsize);
        for i=1:mboxsize
            for j=1:nboxsize
                %              del2eta(i,j)=1/delx^2*(etani(indg(i+1,mboxsize),j)-2*etani(i,j)+etani(indg(i-1,mboxsize),j))...
                %                +1/delx^2*(etani(i,indg(j+1,nboxsize))-2*etani(i,j)+etani(i,indg(j-1,nboxsize)));
                sumeta=etani(indg(i+1,mboxsize),j)+etani(indg(i-1,mboxsize),j)+etani(i,indg(j+1,nboxsize))+etani(i,indg(j-1,nboxsize));
                sumeta2=etani(indg(i+1,mboxsize),indg(j+1,nboxsize))+etani(indg(i-1,mboxsize),indg(j+1,nboxsize))+etani(indg(i+1,mboxsize),indg(j-1,nboxsize))+etani(indg(i-1,mboxsize),indg(j-1,nboxsize));
                del2eta_1(i,j)=2/3/(delx*delx)*(sumeta+0.25*sumeta2-5*etani(i,j));
            end
        end
        del2eta=1*del2eta_1;
        gradeta=-0.5*sqrt(m/2/kappa)*(1-(2*etani-1).^2);
        grad2eta=-sqrt(m/2/kappa)^2*(2*etani-1).*(1-(2*etani-1).^2);
        H=(del2eta-grad2eta)./gradeta;
        H(abs(H)>0.07)=nan;H(etani>0.8)=nan;H(etani<0.2)=nan; % filtering H
        %% ploting
%         gcf;clf;
%         surf(H)
%         view([0 90])
%         shading flat
%         colorbar;axis equal; box on
%         xlabel('X');ylabel('Y'); title([  'Grain number ' num2str(gi)])
%         ax=get(gcf,'children');set(ax(2),'Clim',[-0.03 0.03])
%         pause(0.1)
%         MH(:,:,gni)=H;
        %% sampling based on middle poi2nt of the segment
        BWsegs=~isnan(H);
        Lsegs=bwlabel(BWsegs,8); % 8 neighbor
        prop=regionprops(Lsegs,'Area','BoundingBox');
        Hsamples=[];
        MKSegs=zeros(1,length(prop));
        for segi=1:length(prop)
            % find minimum of H in the segment
            Hseg=H(fix(prop(segi).BoundingBox(2))+1:fix(prop(segi).BoundingBox(2))+fix(prop(segi).BoundingBox(4)), fix(prop(segi).BoundingBox(1))+1:fix(prop(segi).BoundingBox(1))+fix(prop(segi).BoundingBox(3)));
            [Hsegmin]=min(min(Hseg));
            [Hsegmini,Hsegminj]=find(Hseg==Hsegmin);
            selbnd=[fix(0.1*prop(segi).BoundingBox(4)) fix(0.1*prop(segi).BoundingBox(3))]+1; % data selection bounding box is 0.2 of segment bounding box
            try
            Ksegi=Hseg(Hsegmini-selbnd(1):Hsegmini+selbnd(1), Hsegminj-selbnd(2):Hsegminj+selbnd(2));
             nanelems=sum(sum(isnan(Ksegi)));
            Ksegi(isnan(Ksegi))=0;
            Ksegi=sum(sum(Ksegi))/(size(Ksegi,1)*size(Ksegi,2)-nanelems); % average of the Ksegi matrix without considering nan elements
            catch
                Ksegi=Hsegmin; % if the segment is too small then take the min for the whole segment
            end
            if isnan(Ksegi)
                disp('Warning: A too-much-urved interface!!!')
            end
            Hsamples=[Hsamples ones(1,prop(segi).Area)*Ksegi]; % repeats the Ksegi as much as the area of the segment
            MKSegs(segi)=Ksegi;
        end
        
         %%  sampling of H from the whole structure
%         [i,j]=find(~isnan(H)); % converting the H matrix with lots of NaNs to a vercor (Hsamples) for just one grain
%         Hsamples=zeros(1,length(i));
%         for n=1:length(i)
%             Hsamples(n)=H(i(n),j(n));
%         end
        MHsamples{gni}=Hsamples; % keeps all the grains data
    end
    
    %% analysing MHSamples (curvature histograms)
    AllHsamples=[]; % for all the structure conbining all inviditual grains curvarue data to Hsamples
    for ii=1:length(MHsamples)
        AllHsamples=[AllHsamples MHsamples{ii}];
    end
    %     figure
    %     hist (AllHsamples,100)
    Mtn_K{tni}=AllHsamples;
    MHsamples={};
end
clear eta inds
save([ savedir 'curvature_data.mat'])

%%
