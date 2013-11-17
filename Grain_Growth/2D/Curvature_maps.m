clear
savedir='/home/magnetadmin/Documents/Results/2D/Fric2000_m2_k4_Pz20_run2/';
start=8000;
step=8000;
ending=8000;
mboxsize=2000;
nboxsize=2000;
m=2;
kappa=3;
delx=1;
delt=0.1;
pmax=5;
eta=zeros(mboxsize,nboxsize,pmax);
inds=eta;
% tn=40000
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
        H(abs(H)>0.06)=nan;H(etani>0.8)=nan;H(etani<0.2)=nan; % filtering H
        %% ploting
        %         gcf;clf; % figure
        %         surf(H)
        %         view([0 90]); shading flat; colorbar;axis equal; box on
        %         xlabel('X');ylabel('Y'); title([  'Grain number ' num2str(gi)])
        %         ax=get(gcf,'children');set(ax(2),'Clim',[-0.03 0.03])
        %         pause(0.1)
        %         MH(:,:,gni)=H;
        %% sampling based on middle point of the segment
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
% save([ savedir 'curvature_data.mat'])
%%
bigly bigly
%%
clear
savedir='/media/09004e3f-3e8e-409e-b718-5a0ca5495abd/home/cenna/Results/2D/Fric1200s2_m2_k3_phi/10/';
load([ savedir 'curvature_data.mat'])
L=1;
m=2;
kappa=3;
Pz=0.010;
mobility=3/2*L*sqrt(2*kappa/m)
intenergy=1/3*sqrt(2*m*kappa)

%% analysing MHSamples (curvature histograms)
tni=0;
figure
for tn=[start:step:ending]
    tni=tni+1;
    Hsamples=[];
    Hsamples=Mtn_K{tni};
    %     Hsamples(isnan(Hsamples))=[];
    Hsamples((Hsamples)<0)=[]; % additional removing of high curvaure areas
    Hsamples(abs(Hsamples)>0.06)=[];
    %     figure
    %     [N,x]=hist((Hsamples),100);
    Mtn_meanK(tni)=mean(Hsamples);
    Mtn_meanR(tni)=mean(1./Hsamples);

    x=linspace(-10,-2,25); % bins
    [N]=hist (log(abs(Hsamples)),x); % /log(mean(1./Hsamples))
    bar(x,N/length(Hsamples));
    axis([min(x) max(x) 0 0.2])
    title([ 'Time step = ' num2str(tn) ], 'FontSize',16);h=get(gcf,'Children'); set(h(1),'LineWidth',1.5);set(h,'FontSize',14)
    xlabel('ln(Mean Curvarure), ln(K)','FontSize',16); ylabel('Length Fraction','FontSize',16);% pause(1)
    hold on
    plot(log([Mtn_meanK(tni) Mtn_meanK(tni)]), [0 1],'g-' ,'LineWidth',2)
    plot(log([Pz/intenergy Pz/intenergy]), [0 1],'r-' ,'LineWidth',2)
    hold off
    pause(0.05)
    %      mkdir([savedir 'curvarure_dist/'])
    %      print([savedir 'curvarure_dist/' num2str(tni) '.png'],'-dpng','-r200',gcf)
end
%%
Mtn=[start:step:tn];
timevec=Mtn*delt;

figure
plot(timevec,(1./Mtn_meanK).^2,'-o')
% plot(timevec,Mtn_meanR.^2,'rs')
title(savedir)
xlabel('Time')
ylabel('(Curvature Based Equivalent Radius)^2 (1/K)^2')

%% load grain stats
Mtn=[start:step:10000000];
savedir='/media/09004e3f-3e8e-409e-b718-5a0ca5495abd/home/cenna/Results/2D/Fric1200s2_m2_k3_phi/10/';
try
    for tni=1:size(Mtn,2)
        tn=Mtn(tni);
        grains=importdata([savedir 'GrainStat_' num2str(tn) '.txt']);
        grains(grains<5)=[];
        Rgrains=sqrt(grains/pi);
        Rbar(tni)=mean(Rgrains)*delx;
        numgrains(tni)=length(grains);
        biggrain(tni)=grains(1)*delx^2;
    end
catch
end
Mtn=[start:step:tn-step];
timevecG=Mtn*delt;

figure
% plotyy(timevec,(1./Mtn_meanK), timevecG, Rbar)
plot(1./Rbar(1:length(Mtn_meanK)),Mtn_meanK,':.')
ylabel('<K>')
xlabel('1/<R>')
%% comparison between grain stat and curvature (alpha stuff)
timevec=timevec;
alpha=interp1(timevecG, Rbar, timevec).*Mtn_meanK
figure
plot(timevec,alpha,'s')
ylabel( '\alpha (Average Curvature \times Average Radius)')
xlabel(' Time')

%%
beta=Mtn_meanK(end)*intenergy/Pz
alpha/beta

%%
P_z=[0.015 0.030 0.02 0.030 0.010 0.05]
alpha=[0.256 0.239 0.2652 0.2399 0.2908 0.2907]
beta=[0.3347 0.3523 0.3278 0.3457 0.4321 0.3451]

P_z=[0.015 0.030 0.02 0.030  0.05]
alpha=[0.256 0.239 0.2652 0.2399 0.2907]
beta=[0.3347 0.3523 0.3278 0.3457 0.3451]

alphaObeta=alpha./beta

figure
subplot(3,1,1)
plot(P_z, alpha,'o')
ylabel('\alpha')
subplot(3,1,2)
plot(P_z,beta,'rs')
ylabel('\beta')
subplot(3,1,3)
plot(P_z,alpha./beta,'bd')
ylabel('\alpha / \beta')
xlabel('P_z')


%% making a big surf plot from MH
MMH=zeros(mboxsize,nboxsize)*nan;
for i=1:mboxsize
    for j=1:nboxsize
        for p=1:size(MH,3)
            if MH(i,j,p)>0
                MMH(i,j)=MH(i,j,p);
            end
        end
    end
end
figure
surf(phi)
view([0 90])
shading flat
colorbar;axis equal; box on
xlabel('X');ylabel('Y');



%% ------------------------------- Close looks---------------
%% test for a shrinking circle
clear
L=1*1;
m=2;
gamma=1.5*m;
kappa=4;
% geometry settings
p=2;
r=fix(50);
mboxsize=(fix(r)+20); % y axis in pixels
nboxsize=(fix(r)+20); % x axisk
delx=1;      % length unit per pixel
delt=0.03;
% *** Phase Field Procedure *** (so small and simple piece of code!)
eta=zeros(mboxsize,nboxsize,p);
% making initial structure
eta(:,:,1)=circlegrain(mboxsize,nboxsize,nboxsize/2,mboxsize/2,r,'circ');
% eta(:,1:nboxsize/2,1)=1;
eta(:,:,2)=imcomplement(eta(:,:,1));
eta2=eta;
%initialization
for tn=1:300
    for i=1:mboxsize
        for j=1:nboxsize
            %             del2=1/delx^2*(eta(indg(i+1,mboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-1,mboxsize),j,:))...
            %                 +1/delx^2*(eta(i,indg(j+1,nboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-1,nboxsize),:));
            sumeta=eta(indg(i+1,mboxsize),j,:)+eta(indg(i-1,mboxsize),j,:)+eta(i,indg(j+1,nboxsize),:)+eta(i,indg(j-1,nboxsize),:);
            sumeta2=eta(indg(i+1,mboxsize),indg(j+1,nboxsize),:)+eta(indg(i-1,mboxsize),indg(j+1,nboxsize),:)+eta(indg(i+1,mboxsize),indg(j-1,nboxsize),:)+eta(indg(i-1,mboxsize),indg(j-1,nboxsize),:);
            del2=2/3/(delx*delx)*(sumeta+0.25*sumeta2-5*eta(i,j,:));
            sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
            detadtM=m*(-eta(i,j,:)+eta(i,j,:).^3)+2*gamma*sumterm-kappa.*del2;
            detadt=-L.*(detadtM);
            Mdetadt(i,j)=detadt(1);
            eta2(i,j,:)=eta(i,j,:)+delt*detadt;
        end
    end
    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);
    %     drawgrains(phi,nan,nan,tn)
end
%% numerical circle
clear
m=1;
kappa=4;
a=sqrt(m/2/kappa);
r0=20;
delx=1;
syms r
e=1/2*(1-tanh(a*(r-r0)));
% ezplot(e, [0 50])
x=linspace(-35,35,70);
y=x;
[X,Y]=meshgrid(x,y);
e=1/2*(1-tanh(a*((X.^2+Y.^2).^0.5-r0)));
%% calculating curvature
 etani=eta(:,:,1);
% etani=e;
% del2eta=4*del2(etani);
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
H(abs(H)>0.1)=nan;
H(etani>0.9)=nan;
H(etani<0.1)=nan;

%% another way for calculating curvature
etani=eta(:,:,1);
for i=1:mboxsize
    for j=1:nboxsize
        eta_y(i,j)=1/2/delx*(etani(indg(i+1,mboxsize),j)-etani(indg(i-1,mboxsize),j));
        eta_x(i,j)=1/2/delx*(etani(i,indg(j+1,nboxsize))-etani(i,indg(j-1,nboxsize)));
        eta_yy(i,j)=1/delx^2*(etani(indg(i+1,mboxsize),j)-2*etani(i,j)+etani(indg(i-1,mboxsize),j));
        eta_xx(i,j)=1/delx^2*(etani(i,indg(j+1,nboxsize))-2*etani(i,j)+etani(i,indg(j-1,nboxsize)));
        %             eta_yy(i,j)=1/12/delx^2*(-etani(indg(i+2,mboxsize),j)+16*etani(indg(i+1,mboxsize),j)-30*etani(i,j)+16*etani(indg(i-1,mboxsize),j)-etani(indg(i-2,mboxsize),j));
        %             eta_xx(i,j)=1/12/delx^2*(-etani(i,indg(j+2,nboxsize))+16*etani(i,indg(j+1,nboxsize))-30*etani(i,j)+16*etani(i,indg(j-1,nboxsize))-etani(i,indg(j-2,nboxsize)));
        eta_xy(i,j)=1/4/delx^2*(etani(indg(i+1,mboxsize),indg(j+1,nboxsize))-etani(indg(i-1,mboxsize),indg(j+1,nboxsize))-etani(indg(i+1,mboxsize),indg(j-1,nboxsize))+etani(indg(i-1,mboxsize),indg(j-1,nboxsize)));
    end
end
%  for i=1:mboxsize
%         for j=1:nboxsize
%             eta_xy(i,j)=1/4/delx*(eta_y(i,indg(j+1,nboxsize))-eta_y(i,indg(j-1,nboxsize)))+1/4/delx*(eta_x(indg(i+1,mboxsize),j)-eta_x(indg(i-1,mboxsize),j));
%         end
%  end
A=-(eta_x.^2+eta_y.^2);
B=eta_x.^2.*eta_yy+eta_y.^2.*eta_xx-2*eta_x.*eta_y.*eta_xy;
H=-1/1./(sqrt(eta_x.^2+eta_y.^2)).*B./A;
H(abs(H)>0.1)=nan;
H(etani>0.9)=nan;
H(etani<0.1)=nan;
%%
figure
subplot(2,1,1)
plot(eta(:,nboxsize/2),'.-')
hold on
xi=[0:mboxsize-1]*delx; x1=35;
analyticeta=0.5*(1+tanh(sqrt(m/2/kappa)*(x1-xi)));
plot(analyticeta,'r.-')
title(num2str(tn));ylabel('\eta')
hold off
subplot(2,1,2)
plot(H(:,nboxsize/2),'.-'); hold on;
hold off; ylabel('Local Curvature'); xlabel('Position');grid on; axis([0 mboxsize -0.05 -0.03])

%%
figure
surf(H)
view([0 90])
shading flat
colorbar;axis equal; box on
xlabel('X');ylabel('Y')
%
figure
surf(etani)
view([0 90])
shading flat
colorbar;axis equal;
xlabel('X');ylabel('Y')

%% statistic sampling of curvature map

[i,j]=find(~isnan(H));
for n=1:length(i)
    Hsamples(n)=H(i(n),j(n));
end
mean(Hsamples)
figure
hist(Hsamples)
xlabel('Curvature Data')
ylabel('Counut')

%% symbolic derivation of grad2 eta
syms x a
e=0.5*(1+tanh(a*x))
diff(diff(e))

%% analytical circle and its curvature
clear
m=1;
kappa=4;
a=sqrt(m/2/kappa);
r0=20;
syms r
e=1/2*(1-tanh(a*(r-r0)));
% ezplot(e, [0 50])
syms x y
ecirc=subs(e,'r','sqrt(x^2+y^2)')
% ezmesh(ecirc,[-50 50],[-50 50])
eta_x=diff(ecirc,'x');
eta_y=diff(ecirc,'y');
eta_xx=diff(eta_x,'x');
eta_yy=diff(eta_y,'y');
eta_xy=diff(eta_x,'y');
eta_yx=diff(eta_y,'x');


A=-(eta_x^2+eta_y^2);
B=eta_x^2*eta_yy+eta_y^2*eta_xx-2*eta_x.*eta_y.*eta_xy;
H1=simple(-1/(sqrt(eta_x^2+eta_y^2))*B/A)
% ezmesh(H1,[-50 50],[-50 50])
%%
gradeta=-0.5*a*(1-(2*ecirc-1)^2);
grad2eta=-a^2*(2*ecirc-1)*(1-(2*ecirc-1)^2);
del2eta=eta_xx+eta_yy;
H2=simple((del2eta-grad2eta)/gradeta)
% ezmesh(H,[-50 50],[-50 50])
%%
x=linspace(-35,35,70);
y=x;
[X,Y]=meshgrid(x,y);
num_H1=zeros(size(X))*nan;num_H2=num_H1;
for i=1:length(x);
    for j=1:length(y);
        aaa=subs(ecirc,'x',X(i,j));
        num_eta(i,j)=subs(aaa,'y',Y(i,j));
        if num_eta(i,j)>0.1 && num_eta(i,j)<0.9
            aaa=subs(H1,'x',X(i,j));
            num_H1(i,j)=subs(aaa,'y',Y(i,j));
            aaa=subs(H2,'x',X(i,j));
            num_H2(i,j)=subs(aaa,'y',Y(i,j));
        end
    end
end

figure
surf(num_H1)
view([0 90])
shading flat
colorbar;axis equal; box on
xlabel('X');ylabel('Y')


