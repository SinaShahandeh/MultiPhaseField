clear all
savedir='onefield3D_e2';
mkdir(savedir)
% figure;
% phase field parameters
L=[1 1];
alpha=[1 1];
beta=[1 1];
gamma=2;
kappa=[2 2];
epsilon=2;
G=[0 0];
% geometry settings
p=2;
global nboxsize mboxsize lboxsize
global delx
scale=2;
mboxsize=16*scale;%y-direction
nboxsize=20*scale;%x-direction
lboxsize=20*scale;%z-direction
delx=2/scale;
grainD=15*scale;

endtime=100;
timestepn=3000;
delt=endtime/timestepn;


% *** Phase Field Procedure *** (so small and simple piece of code!)
eta=zeros(mboxsize,nboxsize,lboxsize,p);
% making initial structure
eta(:,:,:,1)=spheregrain(mboxsize,nboxsize,lboxsize,nboxsize/2,2.9*mboxsize/4,lboxsize/2,grainD);
eta(:,:,:,2)=imcomplement(eta(:,:,:,1));
% ppf is phase variable representing particles
%particledistro(nboxsize,mboxsize,particles_number,radius)
[ppf,xparticle,yparticle,zparticle]=particledistro3D(nboxsize,mboxsize,lboxsize,1,5*scale);

%savesettings
save(strcat(pwd,'/',savedir,'/','setings.mat'))

eta2=zeros(mboxsize,nboxsize,lboxsize,p); %pre-assignment

for tn=1:10
    for i=1:mboxsize
        for j=1:nboxsize
            for k=1:lboxsize
                del2=1/delx^2*(0.5*(eta(indg(i+1,mboxsize),j,k,:)-2*eta(i,j,k,:)+eta(indg(i-1,mboxsize),j,k,:))...
                    +0.25*(eta(indg(i+2,mboxsize),j,k,:)-2*eta(i,j,k,:)+eta(indg(i-2,mboxsize),j,k,:)))...
                    +1/delx^2*(0.5*(eta(i,indg(j+1,nboxsize),k,:)-2*eta(i,j,k,:)+eta(i,indg(j-1,nboxsize),k,:))...
                    +0.25*(eta(i,indg(j+2,nboxsize),k,:)-2*eta(i,j,k,:)+eta(i,indg(j-2,nboxsize),k,:)))...
                    +1/delx^2*(0.5*(eta(i,j,indg(k+1,lboxsize),:)-2*eta(i,j,k,:)+eta(i,j,indg(k-1,lboxsize),:))...
                    +0.25*(eta(i,j,indg(k+2,lboxsize),:)-2*eta(i,j,k,:)+eta(i,j,indg(k-2,lboxsize),:)));
                sumterm=eta(i,j,k,:)*sum(eta(i,j,k,:).^2)-eta(i,j,k,:).^3;
                detadtM=(-alpha.*reshape(eta(i,j,k,:),1,p)+beta.*reshape(eta(i,j,k,:),1,p).^3-kappa.*reshape(del2,1,p)+...
                    2*epsilon.*reshape(eta(i,j,k,:),1,p)*ppf(i,j,k))+G;
                detadt=-L.*(detadtM+2*gamma*reshape(sumterm,1,p));
                eta2(i,j,k,:)=eta(i,j,k,:)+reshape(delt*detadt,1,1,1,2);
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

                for pind=1:p
                    if eta2(i,j,k,pind)>1
                        eta2(i,j,k,pind)=1;
                    end
                    if eta2(i,j,k,pind)<0
                        eta2(i,j,k,pind)=0;
                    end
                end
            end
        end
    end
    eta=eta2;
    phi=sum(eta(:,:,:,1:p).^2,4);
    % adding ppf to the phi to make particles positions to 1 inorder to
    % make mapping more clear to see it dosen't do anything with the matrix
    phi=phi+ppf;
    drawgrains3D(phi,xparticle,yparticle,zparticle,tn)

    %     [deletax deletay]=gradient(eta,delx,delx);
    %     gradeta=sqrt(deletax.^2+deletay.^2);
    %     E=-alpha/2*eta.^2+beta/4*eta.^4+kappa*gradeta.^2;
    %     ME(tn)=sum(sum(E))/mboxsize/nboxsize;
    %     drawE(rot90(rot90(E)),xparticle,yparticle,tn,eta,ppf)
    savegrains(eta,xparticle,yparticle,zparticle,tn,savedir)

end
%
