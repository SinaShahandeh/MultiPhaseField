% animation

clear
dirstring='/media/disk/sim_res/Ab2'
load(strcat(dirstring,'/','setings','.mat'))
% figure
ni=0;
i=1;
for tn=[1:50 52:2:100 105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000];
    filename=strcat(dirstring,'/',num2str(tn),'.mat');
    load(filename)
    drawgrains(phi,ppf,xparticle,yparticle,tn)
    pause(0.1)
    % M(i)=getframe;

end


%SUBPLOT figures
clear
dirstring='/media/disk/sim_res/onebig4/'
load(strcat(dirstring,'/','setings','.mat'))
figure
ni=0;
sub=100*[1 2 5 8 10 20].^2;
% making Timehistory
Timehistory(1)=0;
for ti=1:length(delt)-1
    Timehistory(ti+1)=Timehistory(ti)+delt(ti);
end
for subi=1:6
    subplot(3,2,subi)
    filename=strcat(dirstring,'/',num2str(sub(subi)),'.mat');
    load(filename)
    drawgrains(phi,ppf,xparticle,yparticle,Timehistory(sub(subi)))
end



%% Size Distribution  %%
clear
savedir='big_grain'
load(strcat(pwd,'/',savedir,'/','setings','.mat'))
% figure
tn=2000
filename=strcat(pwd,'/',savedir,'/',num2str(tn),'.mat');
load(filename)


for i=1:mboxsize
    for j=1:nboxsize
        % calculation of nabla square eta
        del2=1/delx^2*(0.5*(eta(indg(i+1,nboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-1,nboxsize),j,:))...
            +0.25*(eta(indg(i+2,nboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-2,nboxsize),j,:)))...
            +1/delx^2*(0.5*(eta(i,indg(j+1,mboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-1,mboxsize),:))...
            +0.25*(eta(i,indg(j+2,mboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-2,mboxsize),:)));
        % double summation part in the PDE equation. cube term is to
        % substract the etai*eta^2 part from sum and get the j~=i
        % summation
        sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
        detadtM=(-alpha*eta(i,j,:)+beta*eta(i,j,:).^3-kappa*del2+...
            2*epsilon*eta(i,j,:)*ppf(i,j));
        detadt=-L*(detadtM+2*gamma*(sumterm));
        eta2(i,j,:)=eta(i,j,:)+delt(tn)*detadt;
        % for making sure eta is not outside the equilibrium values
        % actually it is unnecessary
        for pind=1:p
            if eta2(i,j,pind)>1
                eta2(i,j,pind)=1;
            end
            if eta2(i,j,pind)<0
                eta2(i,j,pind)=0;
            end
        end

    end
end

[nablaetax,nablaetay]=gradient(eta(:,:,1),delx,delx);
%     nablaetax((phi>0.99))=nan;
%     nablaetay((phi>0.99))=nan;
Mvelocx=Mdetadt./nablaetax;
Mvelocy=Mdetadt./nablaetay;
%         Mveloc=sqrt(Mvelocx.^2+Mvelocy.^2);
Mveloc=Mdetadt./sqrt(nablaetax.^2+nablaetay.^2);


% Image Processing of Phi

%% analysing phi field
clear
dirstring='/media/disk/sim_res/onebig2/'
load(strcat(dirstring,'/','settings','.mat'))
% figure
ni=0;
i=1;
for tn=[105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000];
    filename=strcat(dirstring,'/',num2str(tn),'.mat');
    load(filename)
    phi(logical(ppf))=0;
    [areaG,diamG,perimG]=grainstat(phi);
%     [areaG,diamG,perimG]=grainstat_simple(phi);
    A{i}=areaG;
    D{i}=diamG;
    P{i}=perimG;
    i=i+1;
    tn
end
save(strcat(dirstring,'/','GrainsStatistics','.mat'))


% reading the statistics and analysis
clear
dirstring='/media/disk/sim_res/onebig1/'
load(strcat(dirstring,'/','GrainsStatistics','.mat'))
i=1;
tn=[1:2:100 105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000];
try
    for tni=tn
        Ai=A{i};Di=D{i};Pi=P{i};
        % removing small dots
        [removeindex]=find(Ai<exp(1));
        Ai(removeindex)=[];
        Di(removeindex)=[];
        Pi(removeindex)=[];
        % big grain properties
        [maxA(i),indmaxG]=max(Ai);
        maxD(i)=Di(indmaxG);
        maxP(i)=Pi(indmaxG);
        % removing big grain from statistics of others
        Ai(indmaxG)=[];
        Di(indmaxG)=[];
        Pi(indmaxG)=[];

        An(i)=mean(Ai);
        %   Dn(i)=sqrt(mean((Di).^2));
        Dn(i)=mean(Di);
        Pn(i)=mean(Pi);
        GBArea(i)=sum(Pi)/2;
        i=i+1;
    end
catch
    disp('Warning: simulation was incomplete')
end
% making Timehistory
Timehistory(1)=0;
for ti=1:length(delt)-1
    Timehistory(ti+1)=Timehistory(ti)+delt(ti);
end

figure
plot(Timehistory(tn(1:size(An,2))),An)
hold on
plot(Timehistory(tn(1:size(An,2))),maxA,'r')
title('Area of Abnormal (red) and matrix (blue) Grains,')


plot(tn(1:end),Pn(1:end))
plot(tn(1:end),GBArea(1:end))


%% ----------------------------------------------------------------
%% Distribution Analysis
clear
dirstring='/Drive2/sim_res/Ab2'
load(strcat(dirstring,'/','GrainsStatistics','.mat'))
i=1;
tn=[1:50 52:2:100 105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000];% 20200:200:40000];

GrainsA=A{:};
figure
[n,A]=hist(log(GrainsA),120);
bar(A,n/sum(n))
title('Number density histogram')
xlabel('Ln(Grain Area)')
ylabel('Frequency')

figure
bar(A,n.*exp(A)/sum(n.*exp(A)))
title('Area density histogram')
xlabel('Ln(Grain Area)')
ylabel('Frequency')

%% Extract two size population density Based on number density
% A curve fitting code, fits double Gaussian fucntion to distribution
f_ = figure;
figure(f_);
set(f_,'Units','Pixels','Position',[441 134 680 475]);
legh_ = []; legt_ = {};   % handles and text for legend
xlim_ = [Inf -Inf];       % limits of x axis
ax_ = axes;
set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
set(ax_,'Box','on');
axes(ax_); hold on;
% --- Plot data originally in dataset "y vs. x"
x = A(:);
y = n(:)/sum(n(:));
h_ = line(x,y,'Parent',ax_,'Color',[0.333333 0 0.666667],...
    'LineStyle','none', 'LineWidth',1,...
    'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(x));
xlim_(2) = max(xlim_(2),max(x));
legh_(end+1) = h_;
legt_{end+1} = 'y vs. x';
% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
    xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
    set(ax_,'XLim',xlim_)
else
    set(ax_, 'XLim',[-15.161660148229525902, -3.0549015239848000824]);
end
% --- Create fit "fit 2"
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[-Inf -Inf    0 -Inf -Inf    0]);
ok_ = isfinite(x) & isfinite(y);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [0.051617003725565768768 -8.7591816831589301984 0.33214321585755812727 0.036308748919838901292 -8.1607259923905299814 0.45268087450771804514 ];
set(fo_,'Startpoint',st_);
ft_ = fittype('gauss2');
% Fit this model using new data
cf_ = fit(x(ok_),y(ok_),ft_,fo_);
% Plot this fit
h_ = plot(cf_,'fit',0.95);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'Gaussian Fit';
% Done plotting data and fits.  Now finish up loose ends.
hold off;
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'};
h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
set(h_,'Interpreter','none');
xlabel(ax_,'');               % remove x label
ylabel(ax_,'');               % remove y label
% plotting individual distributions
hold on
coeff=coeffvalues(cf_);
xi=linspace(xlim_(1),xlim_(2),500);
y1=coeff(1)*exp(-((xi-coeff(2))/coeff(3)).^2);
y2=coeff(4)*exp(-((xi-coeff(5))/coeff(6)).^2);
plot(xi,y1,'g-');
plot(xi,y2,'y-.');
grid on

%% reconstructing area frequency density functions
ya1=y1.*sum(n).*exp(xi)/sum(y1.*sum(n).*exp(xi));
ya2=y2.*sum(n).*exp(xi)/sum(y2.*sum(n).*exp(xi));
yatot=ya1+ya2;
figure
plot(xi,yatot,'r')
hold on
plot(xi,ya1,'g-')
plot(xi,ya2,'y-.');
grid on




%% -----------  velocity field:

clear
savedir='Ab1'
load(strcat(pwd,'/',savedir,'/','setings','.mat'))
% figure
ni=0;
i=1;
tns=[1:50 52:2:100 105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000];
tn=400
filename=strcat(pwd,'/',savedir,'/',num2str(tn),'.mat');
load(filename)

eta1=eta;
tn2=tns(find(tns==tn)+1)
filename=strcat(pwd,'/',savedir,'/',num2str(tn2),'.mat');
load(filename)
eta2=eta;

% gradient matrixes
pn=2
[nablaetax,nablaetay]=gradient(eta1(:,:,pn),delx,delx);
% deta/dt field
Mdetadt=(eta2(:,:,pn)-eta1(:,:,pn))/delt(tn);

% Mvelocx=Mdetadt./nablaetax;
% Mvelocy=Mdetadt./nablaetay;
% Mveloc=sqrt(Mvelocx.^2+Mvelocy.^2);
Mveloc=Mdetadt./sqrt(nablaetax.^2+nablaetay.^2);
% Mveloc(Mveloc>3)=0;
% Mveloc(Mveloc<-3)=0;
% Mveloc(Mveloc==nan)=0;
phi(ppf==1)=1;
Mdetadt(ppf==1)=0;
maxM=0.009*max(max(abs(Mdetadt)));
Mveloc(abs(Mdetadt)<maxM)=0;


% plotting velocity field
x=linspace(0,nboxsize*delx,nboxsize);
y=linspace(0,mboxsize*delx,mboxsize);
[X,Y]=meshgrid(x,y);

% h=surf(X,Y,rot90(rot90(MV)));
h=surf(X,Y,((Mveloc)));
shading interp
view([0 90])
% set(h,'clim',[0 1])
colormap jet
% axis equal
% set(get(h,'parent'),'clim',[-1.5 1.5])
colorbar
axis([0 nboxsize*delx 0 mboxsize*delx])



%% For all phase fields
etavel=zeros(nboxsize,mboxsize);
for pn=1:p
    % gradient matrixes
    [nablaetax,nablaetay]=gradient(eta1(:,:,pn),delx,delx);
    % deta/dt field
    Mdetadt=(eta2(:,:,pn)-eta1(:,:,pn))/delt(tn);
    %     Mvelocx=Mdetadt./nablaetax;
    %     Mvelocy=Mdetadt./nablaetay;
    %         Mveloc=sqrt(Mvelocx.^2+Mvelocy.^2);
    Mveloc=Mdetadt./sqrt(nablaetax.^2+nablaetay.^2);
    % Mveloc(Mveloc>3)=0;
    % Mveloc(Mveloc<-3)=0;
    % Mveloc(Mveloc==nan)=0;
    phi(ppf==1)=1;
    Mdetadt(ppf==1)=0;
    maxM=0.002*max(max(abs(Mdetadt)));
    Mveloc(abs(Mdetadt)<maxM)=0;

    etavel=etavel+eta(:,:,pn).*Mveloc;
end


% plotting velocity field
x=linspace(0,nboxsize*delx,nboxsize);
y=linspace(0,mboxsize*delx,mboxsize);
[X,Y]=meshgrid(x,y);

% h=surf(X,Y,rot90(rot90(MV)));
h=surf(X,Y,((etavel)));
shading interp
view([0 90])
% set(h,'clim',[0 1])
colormap jet
axis equal
% set(get(h,'parent'),'clim',[-1.5 1.5])
colorbar
axis([0 nboxsize*delx 0 mboxsize*delx])



