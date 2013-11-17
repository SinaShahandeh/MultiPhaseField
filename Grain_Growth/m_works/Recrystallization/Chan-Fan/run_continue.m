% function runtime_continue
clear
% start after: 
tn=800

% dirstring='/Drive2/sim_res/non_uniform_par1/';
% dirstring='/Drive2/sim_res/Particles1/Dissolve800/';
% dirstring='/Drive2/sim_res/non_uniform_par1/'
dirstring='/media/disk/sim_res/Particles10/Dissolve800/'
load(strcat(dirstring,'setings.mat'))
filename=strcat(dirstring,'D',num2str(tn),'.mat');
load(filename)

savedir=dirstring;
% Changing model parameters
% L=1;
% alpha=1;
% beta=1;
% gamma=1;
% kappa=2;
% epsilon=5;

setings.L=L;
setings.alpha=alpha;
setings.beta=beta;
setings.gamma=gamma;
setings.kappa=kappa;
setings.epsilon=epsilon;


eta2=eta; %pre-assignment
%% changing particles distribution 
% diameter=1;
% particles_fraction=0.00;
% % particles number
% particlesn=particles_fraction*nboxsize*nboxsize/diameter^2
% %% changing particles state
% [ppf,xparticle,yparticle]=particledistro(nboxsize,mboxsize,particlesn,diameter,'normal');

% remove small particles
for pn=1:length(xparticle2)
    ppf(yparticle2(pn),xparticle2(pn))=0;
end


% % changing save time steps
savetimesteps=[1:2:100 105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000];
savedetaile=[50 100:50:500 600:100:3000 3000:200:4000 4500:500:10000 11000:1000:40000];
% timestepn=1460
se=strel('square',3);
while tn<timestepn
    tn=tn+1
    tic
    % findig nodes which are in the grain boundaries and solve
    % differential equation only for that points.
    [yii,xjj]=find(...
        imerode((phi>0.999),se)==0);
    % space discretization loop
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

    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);
    % adding ppf to the phi to make particles positions to 1 inorder to
    % make mapping more clear to see. Because the range of phi changes
    % from 0.5 to 1 and having zero element makes range broader and
    % lower contrast. It dosen't do anything with the eta matrix
    %     phi=phi+ppf;
    % draw gray scale image of structure. with particles in red dots
%   drawgrains(phi,ppf,xparticle,yparticle,tn)
  drawgrains(phi,0,0,0,tn)
    %% saving the structure
    % desired saving time steps
    if ~isempty(find(savetimesteps==tn)) && saveresults==1
        filename=strcat(savedir,num2str(tn),'.mat');
        save(filename,'phi','tn') % one can add eta, for having all phase parameters
    end
    if ~isempty(find(savedetaile==tn)) && saveresults==1
        filename=strcat(savedir,'D',num2str(tn),'.mat');
        save(filename,'phi','tn','eta','Timehistory') % one can add eta, for having all phase parameters
    end
    %     savegrains(phi,xparticle,yparticle,tn,savedir)
   % display speed of this step calulation. Using
    toc
end

break
%% ------------- STEP 2 of calculation -------------------

L=1;
alpha=1;
beta=1;
gamma=1;
kappa=2;
epsilon=5;
eta2=eta; %pre-assignment

diameter=1;
particles_fraction=0.00;
% particles number
particlesn=particles_fraction*nboxsize*nboxsize/diameter^2
%% changing particles state
[ppf,xparticle,yparticle]=particledistro(nboxsize,mboxsize,particlesn,diameter,'normal');

% % changing save time steps
% savetimesteps=[1:50 52:2:100 105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000];
% savedetaile=[50 100 500 1000 2000 3000 4000 5000 7000 10000 15000 20000 25000 30000 35000];
timestepn=40000
se=strel('square',3);
while tn<timestepn
    tn=tn+1
    tic
    % findig nodes which are in the grain boundaries and solve
    % differential equation only for that points.
    [yii,xjj]=find(...
        imerode((phi>0.999),se)==0);
    % space discretization loop
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

    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);
    % adding ppf to the phi to make particles positions to 1 inorder to
    % make mapping more clear to see. Because the range of phi changes
    % from 0.5 to 1 and having zero element makes range broader and
    % lower contrast. It dosen't do anything with the eta matrix
    %     phi=phi+ppf;
    % draw gray scale image of structure. with particles in red dots
  % drawgrains(phi,ppf,xparticle,yparticle,tn)
    %% saving the structure
    % desired saving time steps
    if ~isempty(find(savetimesteps==tn)) && saveresults==1
        filename=strcat(savedir,num2str(tn),'.mat');
        save(filename,'phi','tn') % one can add eta, for having all phase parameters
    end
    if ~isempty(find(savedetaile==tn)) && saveresults==1
        filename=strcat(savedir,'D',num2str(tn),'.mat');
        save(filename,'phi','tn','eta','Timehistory') % one can add eta, for having all phase parameters
    end
    %     savegrains(phi,xparticle,yparticle,tn,savedir)
   % display speed of this step calulation. Using
    toc
end



