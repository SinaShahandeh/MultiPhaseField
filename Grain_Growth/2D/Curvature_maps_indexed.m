%% curvature maps based on indexed matrix

%% Reading Inds
clear
start=2000;
step=2000;
ending=98000;
tni=0;
m=1;
kappa=4;
a=sqrt(m/2/kappa);
delx=1;
delt=0.05;

tn=20000;
%for tn=[10:10:20000000]
savedir='/home/cenna/Results/2D/Fric250s2_m1_k4/0/'
Inds=importdata([savedir '/Inds_' num2str(tn) '.txt']);
figure
surf(Inds);
view ([0 90])
shading flat
axis off
box off
title(['Time= ' num2str(tn)])
set(gca,'DataAspectratio',[1 1 1])
% end

GrainStat=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
BiggestGrainD=2*sqrt(max(GrainStat)/pi); % to repeat the the system at the boundaries to include the largest grain if it falls on the boundary
BigGrains=find(GrainStat>10)';
gni=0; % counter for grains
for gi=BigGrains
    gni=gni+1;
    mboxsize=size(Inds,1);nboxsize=size(Inds,2);
    indsgi=zeros(mboxsize,nboxsize);
    for i=1:mboxsize
        for j=1:nboxsize
            if (Inds(i,j)==gi)
                indsgi(i,j)=Inds(i,j); % separating grain index gi
            end
        end
    end
    [in,jn]=find(indsgi>0);
    jmin=min(jn);jmax=max(jn);imin=min(in);imax=max(in);
    % Transfering the domain to make the bounudary grains at the centre
    if jmax+BiggestGrainD>nboxsize
        indsgi=[indsgi(:,nboxsize/2:end)  indsgi(:,1:nboxsize/2-1)];
    end
    if jmin-BiggestGrainD<1
        indsgi=[indsgi(:,nboxsize/2:end)  indsgi(:,1:nboxsize/2-1)];
    end
    if imax+BiggestGrainD>mboxsize
        indsgi=[indsgi(mboxsize/2:end,:); indsgi(1:mboxsize/2-1,:)];
    end
    if imin-BiggestGrainD<1
        indsgi=[indsgi(mboxsize/2:end,:); indsgi(1:mboxsize/2-1,:)];
    end
    if (imax==mboxsize) && (imin==1)
        indsgi=[indsgi(mboxsize/2:end,:); indsgi(1:mboxsize/2-1,:)];
    end
    if (jmin==1) && (jmax==nboxsize)
        indsgi=[indsgi(:,nboxsize/2:end)  indsgi(:,1:nboxsize/2-1)];
    end
    % cropping the grain from the whole domain
    [in,jn]=find(indsgi>0);
    jmin=min(jn);jmax=max(jn);imin=min(in);imax=max(in);
    indsgi=indsgi(imin:imax,jmin:jmax);
    mboxsize=size(indsgi,1); nboxsize=size(indsgi,2);
    % creating x{gni} and y{gni} cells that contain all the preferal pixels of
    % grain gi
    pc=0; % point counter
    x=[];y=[];
    for i=1:mboxsize
        rowi=indsgi(i,:); % selecting each row and find which point has gi in it and then find first and last element of that row
        [ii]=find(rowi==gi);
        if ~isempty(ii)
            pc=pc+1;
            x(pc)=ii(1);y(pc)=i;
            pc=pc+1;
            x(pc)=ii(end);y(pc)=i;
        end
    end
    for j=1:nboxsize % do the same thing for y direction
        colj=indsgi(:,j); % selecting each row and find which point has gi in it and then find first and last element of that row
        [jj]=find(colj==gi);
        if ~isempty(jj)
            pc=pc+1;
            x(pc)=j;y(pc)=jj(1);
            pc=pc+1;
            x(pc)=j;y(pc)=jj(end);
        end
    end
    x=x*delx;y=y*delx;
    xcent=mean(x);ycent=mean(y);
    [teta,ro]=cart2pol(x-xcent,y-ycent);
    [teta,IX]=sort(teta);
    ro=ro(IX);
    repindex=[diff(teta)==0];
    teta(repindex)=[];
    ro(repindex)=[];
    polar_sp=csaps(teta,ro,0.9); % r as a function of teta
       
    % curvature calculation
    tetai=linspace(-pi,pi,720);
    ri=fnval(polar_sp,tetai);
    drdteta=fnval(fnder(polar_sp),tetai);
    d2rdteta2=fnval(fnder(polar_sp,2),tetai);
    Kgni=(2*drdteta.^2-d2rdteta2.*ri+ri.^2)./(ri.^2+drdteta.^2).^(3/2);
    removeKind=Kgni>0.06;
    ri(removeKind)=nan;
    tetai(removeKind)=nan;
    Kgni(removeKind)=nan;
    
    % see curvature in a polar plot
    roi=fnval(polar_sp,tetai);
    figure;subplot(1,2,1)
    polar(teta,ro,'.'); hold on
    polar(tetai,ri) % see fitted spline
    subplot(1,2,2)
    polar(tetai,Kgni,'r')
    
end
