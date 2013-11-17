clear all
savedir='onefield';
mkdir(savedir)
% figure;
clf
% phase field parameters
L=[1];
alpha=[1];
beta=[1];
gamma=2;
kappa=[2];
epsilon=5;

% geometry settings
p=1;
global nboxsize mboxsize
global delx
scale=2;
mboxsize=20*scale;
nboxsize=20*scale;
delx=2/scale;
grainD=15*scale;

endtime=100;
timestepn=4000;
delt=0.10;

%savesettings
save(strcat(pwd,'/',savedir,'/','setings.mat'))

% *** Phase Field Procedure *** (so small and simple piece of code!)
eta=zeros(mboxsize,nboxsize,p)-1;
% making initial structure
eta(:,:,1)=(circlegrain(mboxsize,nboxsize,nboxsize/2,3*mboxsize/4,grainD));
% ppf is phase variable representing particles
%particledistro(nboxsize,mboxsize,particles_number,radius)
[ppf,xparticle,yparticle]=particledistro(nboxsize,mboxsize,1,2*scale);
ppf2=imcomplement(ppf);
% eta(:,5,1)=1;
% eta(:,:,2)=1;
% eta(:,5,2)=0;
% eta=rand(gridn,gridn,p)*0.001;%-0.001;
eta2=zeros(mboxsize,nboxsize,p); %pre-assignment

for tn=1:timestepn
    for i=1:mboxsize
        for j=1:nboxsize
            del2=1/delx^2*(0.5*(eta(indg(i+1,mboxsize),j)-2*eta(i,j)+eta(indg(i-1,mboxsize),j))...
                +0.25*(eta(indg(i+2,mboxsize),j)-2*eta(i,j)+eta(indg(i-2,mboxsize),j)))...
                +1/delx^2*(0.5*(eta(i,indg(j+1,nboxsize))-2*eta(i,j)+eta(i,indg(j-1,nboxsize)))...
                +0.25*(eta(i,indg(j+2,nboxsize))-2*eta(i,j)+eta(i,indg(j-2,nboxsize))));
            sumterm=eta(i,j)*sum(eta(i,j).^2)-eta(i,j).^3;
            detadtM=(-alpha.*eta(i,j)+beta.*eta(i,j).^3-kappa.*del2+...
                2*epsilon.*eta(i,j)*ppf(i,j));
            detadt=-L.*(detadtM+2*gamma*sumterm)*ppf2(i,j);
            eta2(i,j)=eta(i,j)+delt*detadt;
            %             elovution of particle
            %             del2=1/delx^2*(0.5*(ppf(indg(i+1,gridn),j)-2*ppf(i,j)+ppf(indg(i-1,gridn),j))...
            %                 +0.25*(ppf(indg(i+2,gridn),j)-2*ppf(i,j)+ppf(indg(i-2,gridn),j)))...
            %                 +1/delx^2*(0.5*(ppf(i,indg(j+1,gridn))-2*ppf(i,j)+ppf(i,indg(j-1,gridn)))...
            %                 +0.25*(ppf(i,indg(j+2,gridn))-2*ppf(i,j)+ppf(i,indg(j-2,gridn))));
            %             sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
            %             detadtM=(-alpha.*reshape(eta(i,j,:),1,p)+beta.*reshape(eta(i,j,:),1,p).^3-kappa.*reshape(del2,1,p)+...
            %                 2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j))+G;
            %             detadt=-L.*(detadtM+2*gamma*reshape(sumterm,1,p));
            %             eta2(i,j,:)=eta(i,j,:)+reshape(delt*detadt,1,1,2);

            if eta2(i,j)>1
                eta2(i,j)=1;
            end
            if eta2(i,j)<-1
                eta2(i,j)=-1;
            end
        end
    end

    eta=eta2;
    %{
    phi=eta.^2;
    % adding ppf to the phi to make particles positions to 1 inorder to
    % make mapping more clear to see it dosen't do anything with the matrix
    phi=phi+ppf;
    draw(phi,xparticle,yparticle,tn,eta,ppf)
    %}
    [deletax deletay]=gradient(eta,delx,delx);
    gradeta=sqrt(deletax.^2+deletay.^2);
    E=-alpha/2*eta.^2+beta/4*eta.^4+kappa*gradeta.^2;
    ME(tn)=sum(sum(E))/mboxsize/nboxsize;
       drawE(rot90(rot90(E)),xparticle,yparticle,tn,eta,ppf)
%     savegrains(eta,ppf,E,xparticle,yparticle,tn,savedir)
end
%

