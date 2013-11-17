%%% This is Multi Phase Field simulation for normal grain growth in 2D [1,2].
%%% The particle pinning based on Molean's functional [3] is also included.
%
% References:
% 1. L-Q Chen, W. Yang, Phys. Rev. B, 1994, 50; 15752.
% 2. D. Fan, L. Q. Chen, Acta Mater. 1997, 45; 611.
% 3. N. Moelans, B. Blanpain, P. Wollants, Acta Mater. 2005, 54;1175.
%
% function realtime_particles()
clear
savedir='/media/disk-3/sim_res/test/';

load(strcat(savedir,'settings.mat'))
load([savedir 'initeta.mat'],'eta')

processn=4; % number of processors
processi=2; % processor ID number


%% start working on its own task
eta=eta(:,(processi-1)*mboxsize/processn+1:processi*mboxsize/processn,:);
eta2=zeros(nboxsize,mboxsize/processn,p); %pre-assignment

for tn=1:1
    tic
    for i=1:mboxsize/processn
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
    eta=eta2;
    phii(:,:,processi)=sum(eta(:,:,1:p).^2,3);
end

%% making phi field from all processors
for proi=2:processn
    load([savedir 'phi' num2str(proi) '.mat'])

end

% strel element for the kernel of imerode
se=strel('square',3);
% time discretization loop
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
   drawgrains(phi,ppf,xparticle,yparticle,tn)
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

