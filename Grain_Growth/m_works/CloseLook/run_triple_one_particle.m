
%%This is for obtaining equilibrium profiles at triple junctions
function [ETriple,EGB,E,Ep]=run_triple_one_particle(gamma,filenum)
savedir='/media/disk/sim_res/triple_variables02/';
mkdir(savedir)
% figure
% phase field parameters
L=5;
alpha=1;
beta=1;
% gamma=1;
kappa=1;
epsilon=5;
DelG=[0 0 0];
% setings structure
settings.L=L;
settings.alpha=alpha;
settings.beta=beta;
settings.gamma=gamma;
settings.kappa=kappa;
settings.epsilon=epsilon;
settings.DelG=DelG;
settings.accuracy='low';
% phase fields geometry settings
p=3;
global nboxsize mboxsize
global delx delt scale
scale=3;
mboxsize=30*scale; % y axis in pixels
nboxsize=28*scale; % x axis
delx=2/scale;      % length unit per pixel
grainD=nboxsize+1;  % in pixels

% Particle geometry
pr=4;  % size in length unit
xparticle=fix(nboxsize/2); % position on the grid
yparticle=fix(mboxsize/2*10/10);

endtime=50;
timestepn=200;
timestepn2=50;
delt=endtime/timestepn;
delt=0.02
% *** Phase Field Procedure *** (so small and simple piece of code!)
%% initial Structure
eta=zeros(mboxsize,nboxsize,p);
%making initial structure
eta(:,:,1)=circlegrain(mboxsize,nboxsize,nboxsize,6*mboxsize/10,grainD,'dome');
eta(:,:,2)=circlegrain(mboxsize,nboxsize,0,6*mboxsize/10,grainD,'dome');
eta(:,:,3)=imcomplement(eta(:,:,1)+eta(:,:,2));


% % triple junction equilibriation!
% eta(1:mboxsize/3,:,1)=ones(mboxsize/3,nboxsize);
% eta(mboxsize/3+1:mboxsize,1:nboxsize/2,2)=ones(mboxsize-mboxsize/3,nboxsize/2);
% eta(mboxsize/3+1:mboxsize,nboxsize/2+1:nboxsize,3)=ones(mboxsize-mboxsize/3,nboxsize/2);
% 


% ppf is phase variable representing particles
%particledistro(nboxsize,mboxsize,particles_number,radius)
[ppf]=particledistro(nboxsize,mboxsize,0,pr/delx*2,xparticle,yparticle);

eta2=zeros(mboxsize,nboxsize,p); %pre-assignment

%savesettings
save(strcat(savedir,'/','settings.mat'))

%initialization
for tn=1:5
    for i=1:mboxsize
        for j=1:nboxsize
            del2=1/delx^2*(0.5*(eta(indgsym(i+1,mboxsize),j,:)-2*eta(i,j,:)+eta(indgsym(i-1,mboxsize),j,:))...
                +0.25*(eta(indgsym(i+2,mboxsize),j,:)-2*eta(i,j,:)+eta(indgsym(i-2,mboxsize),j,:)))...
                +1/delx^2*(0.5*(eta(i,indgsym(j+1,nboxsize),:)-2*eta(i,j,:)+eta(i,indgsym(j-1,nboxsize),:))...
                +0.25*(eta(i,indgsym(j+2,nboxsize),:)-2*eta(i,j,:)+eta(i,indgsym(j-2,nboxsize),:)));
            sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
            detadtM=(-alpha.*reshape(eta(i,j,:),1,p)+beta.*reshape(eta(i,j,:),1,p).^3-kappa.*reshape(del2,1,p)+...
                2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j));
            detadt=-L.*(detadtM+2*gamma(1)*reshape(sumterm,1,p));
            eta2(i,j,:)=eta(i,j,:)+reshape(delt*detadt,1,1,p);
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
    phi=sum(eta(:,:,1:p).^2,3);
    phi=phi+ppf;
end
se=strel('square',1);
timevec=0;
tn=0;
while tn<timestepn
    Mdetadt=zeros(mboxsize,nboxsize);
    tn=tn+1;
    [yii,xjj]=find(...
        imerode((phi>0.99999),se)==0);
    for ii=1:length(xjj)
        i=yii(ii);j=xjj(ii);
        del2=1/delx^2*(0.5*(eta(indgsym(i+1,mboxsize),j,:)-2*eta(i,j,:)+eta(indgsym(i-1,mboxsize),j,:))...
                +0.25*(eta(indgsym(i+2,mboxsize),j,:)-2*eta(i,j,:)+eta(indgsym(i-2,mboxsize),j,:)))...
                +1/delx^2*(0.5*(eta(i,indgsym(j+1,nboxsize),:)-2*eta(i,j,:)+eta(i,indgsym(j-1,nboxsize),:))...
                +0.25*(eta(i,indgsym(j+2,nboxsize),:)-2*eta(i,j,:)+eta(i,indgsym(j-2,nboxsize),:)));
        sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
        detadtM=(-alpha.*reshape(eta(i,j,:),1,p)+beta.*reshape(eta(i,j,:),1,p).^3-kappa.*reshape(del2,1,p)+...
            2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j)); 
%             +2*ppf(i,j)*[eta(i,j,2)*eta(i,j,3) eta(i,j,1)*eta(i,j,3) eta(i,j,1)*eta(i,j,2)];
        
        detadt=-L.*(detadtM+2*gamma(1)*reshape(sumterm,1,p));
        eta2(i,j,:)=eta(i,j,:)+reshape(delt*detadt,1,1,p);

        for pind=1:p
            if eta2(i,j,pind)>1
                eta2(i,j,pind)=1;
            end
            if eta2(i,j,pind)<0
                eta2(i,j,pind)=0;
            end
        end
    end
    
    Mdetadt=(eta2(:,:,1)-eta(:,:,1))/delt;
    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);

    timevec(tn+1)=timevec(tn)+delt;
    
    %%% Visulaizations

    %% simple structure
    %   subplot(2,3,1)
    %
%     drawgrains(phi,xparticle,yparticle,tn)
%     imshow(eta)
    %     contourgrains(eta,xparticle,yparticle,tn,ppf)
    %% Energy plots
        [ME,E]=calculateE(eta,ppf,mboxsize,nboxsize,delx,settings);
    % VecE(tn)=E;
    %     drawE(rot90(rot90(ME)),xparticle,yparticle,tn,eta,ppf)

    %% Draw Energy Field
%     h=subplot(2,3,3);
    drawvelocity(ME)
    % %     set(h,'clim',[-0.2 0])
    title(['Energy density, timestep= ' num2str(tn)])
pause(0.01)

end
%

%% 

save([savedir num2str(filenum) '_without_particle.mat'])
%% meausing maximum Energy Density
ind=findtripple(phi);
ETriple=ME(fix(ind(2))+1,fix(ind(1))+1);
EGB=max(ME(end,:));

%% add the particle at triple junction

[ppf]=particledistro(nboxsize,mboxsize,1,pr/delx*2,fix(ind(1))+1,fix(ind(2))+1);
delt=delt/2
while tn<timestepn+timestepn2
    Mdetadt=zeros(mboxsize,nboxsize);
    tn=tn+1;
    [yii,xjj]=find(...
        imerode((phi>0.99999),se)==0);
    for ii=1:length(xjj)
        i=yii(ii);j=xjj(ii);
        del2=1/delx^2*(0.5*(eta(indgsym(i+1,mboxsize),j,:)-2*eta(i,j,:)+eta(indgsym(i-1,mboxsize),j,:))...
                +0.25*(eta(indgsym(i+2,mboxsize),j,:)-2*eta(i,j,:)+eta(indgsym(i-2,mboxsize),j,:)))...
                +1/delx^2*(0.5*(eta(i,indgsym(j+1,nboxsize),:)-2*eta(i,j,:)+eta(i,indgsym(j-1,nboxsize),:))...
                +0.25*(eta(i,indgsym(j+2,nboxsize),:)-2*eta(i,j,:)+eta(i,indgsym(j-2,nboxsize),:)));
        sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
        detadtM=(-alpha.*reshape(eta(i,j,:),1,p)+beta.*reshape(eta(i,j,:),1,p).^3-kappa.*reshape(del2,1,p)+...
            2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j)); 
%             +2*ppf(i,j)*[eta(i,j,2)*eta(i,j,3) eta(i,j,1)*eta(i,j,3) eta(i,j,1)*eta(i,j,2)];
        
        detadt=-L.*(detadtM+2*gamma(1)*reshape(sumterm,1,p));
        eta2(i,j,:)=eta(i,j,:)+reshape(delt*detadt,1,1,p);

        for pind=1:p
            if eta2(i,j,pind)>1
                eta2(i,j,pind)=1;
            end
            if eta2(i,j,pind)<0
                eta2(i,j,pind)=0;
            end
        end
    end
    
    Mdetadt=(eta2(:,:,1)-eta(:,:,1))/delt;
    eta=eta2;
    phi=sum(eta(:,:,1:p).^2,3);

    timevec(tn+1)=timevec(tn)+delt;
    
    %%% Visulaizations

    %% simple structure
    %   subplot(2,3,1)
    %
%     drawgrains(phi,xparticle,yparticle,tn)
%     imshow(eta)
    %     contourgrains(eta,xparticle,yparticle,tn,ppf)
    %% Energy plots
        [ME,Ep]=calculateE(eta,ppf,mboxsize,nboxsize,delx,settings);
    % VecE(tn)=E;
    %     drawE(rot90(rot90(ME)),xparticle,yparticle,tn,eta,ppf)

    %% Draw Energy Field
%     h=subplot(2,3,3);
    drawvelocity(ME)
    % %     set(h,'clim',[-0.2 0])
    title(['Energy density, timestep= ' num2str(tn)])
    pause(0.01)
end
%

save([savedir num2str(filenum) '_with_particle.mat'])



