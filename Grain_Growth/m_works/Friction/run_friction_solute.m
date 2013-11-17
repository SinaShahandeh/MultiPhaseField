
% function run_friction_solute(savedir,Pz)
 savedir='/media/disk/sim_res/run_friction_solute_01/';
% choose whether save or not
saveresults=1;
% time steps that save command is executed
savetimesteps=[1:5:100 105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000];
savedetaile=[50 100:50:500 600:100:3000 3000:200:4000 4500:500:10000 11000:1000:40000];
mkdir(savedir)
% figure;
% clf
% phase field parameters
L=1;
alpha=1 ;
beta=1 ;
gamma=1.5;
kappa=2;
epsilon=5;

% solte drag parameters
a=0.2;
b=4;

setings.L=L;
settings.alpha=alpha;
settings.beta=beta;
settings.gamma=gamma;
settings.kappa=kappa;
settings.epsilon=epsilon;


% geometry settings
p=36; % phase field numbers
scale=1;
mboxsize=300*scale; % x axis in pixels
nboxsize=300*scale; % y axis
delx=2/scale;      % length unit per pixel

% time discritization. The maximum delta t should not be very larger and
% depends on the value of the delx to maintain the stability of
% calculation
delt=[linspace(0.25,0.25,100) linspace(0.25,0.3,900) linspace(0.3,0.3,3000)...
    linspace(0.3,0.3,6000) linspace(0.3,0.4,10000) linspace(0.4,0.4,20000)];

timestepn=size(delt,2);
% number of nucleas at the beginning of simulation
nuclein= mboxsize*nboxsize/20; % ~5 percent of grid points are nuclei
% particles distribution specification
% diameter=1;
% particles_fraction=0.000;
% particles number
% particlesn=particles_fraction*nboxsize*nboxsize/diameter^2;

%% for making non-uniform distribution of the particles use particledistroN
% ratio of banded phase fraction to uniform distribution outside band
%ratio=50;

% *** Phase Field Procedure *** (so small and simple piece of code!)
eta=zeros(nboxsize,mboxsize,p);
% putting initial nucleas 1:total number of nuclea (so nucleation density
% can be set)
% One can use non-uniform distribution of nuclei to have bimodal structure
for nn=1:nuclein;
    ii=fix(nboxsize*rand(1,1))+1;jj=fix(mboxsize*rand(1,1))+1;
    eta(ii,jj,fix(p*rand(1,1))+1)=1;
end
% non-uniform particles distribution in a banded structure

% save settings and initial conditions
save(strcat(savedir,'settings.mat'))
% eta 2 holds phase field parameters in a time step i+1
eta2=zeros(nboxsize,mboxsize,p); %pre-assignment
Timehistory=zeros(1,timestepn+1); % pre-assignment timehistory
for tn=1:60
    tic
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
            detadtM=(-alpha*eta(i,j,:)+beta*eta(i,j,:).^3-kappa*del2);
            detadt=-L*(detadtM+2*gamma*(sumterm));
            eta2(i,j,:)=eta(i,j,:)+delt(tn)*detadt;
            % for making sure eta is not outside the equilibrium values
            % actually it is unnecessary
        end
    end
    Mdetadt=eta2-eta;
    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);
    % for making particle positions 1
    %     phi=phi+ppf;
 %  drawgrains(phi,nan,nan,nan,tn)
 %% saving the structure
    % desired saving time steps
    if ~isempty(find(savetimesteps==tn)) && saveresults==1
        filename=strcat(savedir,num2str(tn),'.mat')
        save(filename,'phi','tn') % one can add eta, for having all phase parameters
    end
    if ~isempty(find(savedetaile==tn)) && saveresults==1
        filename=strcat(savedir,'D',num2str(tn),'.mat');
        save(filename,'phi','tn','eta','Timehistory') % one can add eta, for having all phase parameters
        imwrite(imadjust(phi),strcat(savedir,num2str(tn),'.png'),'png')
    end
    %     savegrains(phi,xparticle,yparticle,tn,savedir)
    % storing time history vector for the case where delt varies
    Timehistory(tn+1)=Timehistory(tn)+delt(tn);
    % display speed of this step calulation. Using
    toc
    pause(0.1)
end
% strel element for the kernel of imerode
se=strel('square',3);
% time discretization loop
while tn<timestepn
    tn=tn+1;
    tic
    % findig nodes which are in the grain boundaries and solve
    % differential equation only for that points.
    [yii,xjj]=find(...
        imerode((phi>0.999),se)==0);
    % space discretization loop
    % Velocity map
    %         Mdetadt=(eta2(:,:,1)-eta(:,:,1))/delt;
for pi=1:p
    [nablaetax,nablaetay]=gradient(eta(:,:,pi),delx,delx);
    Mveloc=Mdetadt(:,:,pi)./sqrt(nablaetax.^2+nablaetay.^2);
    % remove wrong data from Mveloc
    [maxdetadt,indimaxdetadt]=max(abs(Mdetadt(:,:,pi)));
    [maxdetadt,indjmaxdetadt]=max(abs(maxdetadt));
    indimaxdetadt=indimaxdetadt(indjmaxdetadt); %indeceis of maximum in detadt where we know maximum of velocity is calculated correctly
    Mveloc(abs(Mveloc)>abs(Mveloc(indimaxdetadt,indjmaxdetadt)))=nan;
    Mveloc=abs(Mveloc);
    Pz(:,:,pi)=a*Mveloc./(1+b^2*Mveloc.^2);
end
    Pz(isnan(Pz))=0;

    for ii=1:length(xjj)
        i=yii(ii);j=xjj(ii);
        % calculation of nabla square eta
        del2=1/delx^2*(0.5*(eta(indg(i+1,nboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-1,nboxsize),j,:))...
            +0.25*(eta(indg(i+2,nboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-2,nboxsize),j,:)))...
            +1/delx^2*(0.5*(eta(i,indg(j+1,mboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-1,mboxsize),:))...
            +0.25*(eta(i,indg(j+2,mboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-2,mboxsize),:)));
        % double summation part in the PDE equation. cube term is to
        % substract the etai*eta^2 part from sum and get the j~=i
        % summation
        sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
        detadtM=-alpha*eta(i,j,:)+beta*eta(i,j,:).^3-kappa*del2+...
            +2*gamma*(sumterm);%+2*epsilon*eta(i,j,:)*ppf(i,j);
        detadt=-L*(Pz(i,j,:)<abs(detadtM)).*...
            (detadtM+6*eta(i,j,:).*(1-eta(i,j,:)).*sign(detadtM).*Pz(i,j,:)/2);
        eta2(i,j,:)=eta(i,j,:)+delt(tn)*detadt;
        % for making sure eta is not outside the equilibrium values
        % actually it is unnecessary
    end
    eta(eta>1)=1;
    eta(eta<0)=0;
    Mdetadt=eta2-eta;
    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);
    % draw gray scale image of structure. with particles in red dots
 %  drawgrains(phi,nan,nan,nan,tn)
    %% saving the structure
    % desired saving time steps
    if ~isempty(find(savetimesteps==tn)) && saveresults==1
        filename=strcat(savedir,num2str(tn),'.mat')
        save(filename,'phi','tn') % one can add eta, for having all phase parameters
    end
    if ~isempty(find(savedetaile==tn)) && saveresults==1
        filename=strcat(savedir,'D',num2str(tn),'.mat');
        save(filename,'phi','tn','eta','Timehistory') % one can add eta, for having all phase parameters
         imwrite(imadjust(phi),strcat(savedir,num2str(tn),'.png'),'png')
    end
    %     savegrains(phi,xparticle,yparticle,tn,savedir)
    % storing time history vector for the case where delt varies
    Timehistory(tn+1)=Timehistory(tn)+delt(tn);
    % display speed of this step calulation. Using
    toc
    pause(0.1)
end
filename=strcat(savedir,'final.mat');
save(filename)
%

%% ----------------------------------------------------------
%Dependent functions

