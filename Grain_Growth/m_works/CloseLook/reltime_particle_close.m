clear
% savedir='set1'
% mkdir(savedir)
% figure;
clf
% phase field parameters
L=1;
alpha=1;
beta=1;
gamma=1;
kappa=2;
epsilon=5;
% geometry settings
p=2;
gridn=40;
delx=2;
delt=0.25;
timestepn=2000;

%savesettings
% save(strcat(pwd,'/',savedir,'/','setings.mat'))

% *** Phase Field Procedure *** (so small and simple piece of code!)
eta=zeros(gridn,gridn,p);
% making initial structure
grainD=20;
eta(:,:,1)=circlegrain(gridn,gridn,gridn/2,gridn/2,grainD);
eta(:,:,2)=imcomplement(eta(:,:,1));
% ppf is phase variable representing particles
%particledistro(nboxsize,mboxsize,particles_number,radius)
[ppf,xparticle,yparticle]=particledistro(gridn,gridn,1,1);
% eta(:,5,1)=1;
% eta(:,:,2)=1;
% eta(:,5,2)=0;
% eta=rand(gridn,gridn,p)*0.001;%-0.001;
eta2=zeros(gridn,gridn,p); %pre-assignment

for tn=1:timestepn
    for i=1:gridn
        for j=1:gridn
            del2=1/delx^2*(0.5*(eta(indg(i+1,gridn),j,:)-2*eta(i,j,:)+eta(indg(i-1,gridn),j,:))...
                +0.25*(eta(indg(i+2,gridn),j,:)-2*eta(i,j,:)+eta(indg(i-2,gridn),j,:)))...
                +1/delx^2*(0.5*(eta(i,indg(j+1,gridn),:)-2*eta(i,j,:)+eta(i,indg(j-1,gridn),:))...
                +0.25*(eta(i,indg(j+2,gridn),:)-2*eta(i,j,:)+eta(i,indg(j-2,gridn),:)));
            sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
            detadtM=(-alpha*eta(i,j,:)+beta*eta(i,j,:).^3-kappa*del2+...
                2*epsilon*eta(i,j,:)*ppf(i,j).^2);
            detadt=-L*(detadtM+2*gamma*(sumterm));
            eta2(i,j,:)=eta(i,j,:)+delt*detadt;
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
    % adding ppf to the phi to make particles positions to 1 inorder to
    % make mapping more clear to see it dosen't do anything with the matrix
    phi=phi+ppf;
   drawgrains(phi,xparticle,yparticle,tn)
%     savegrains(phi,xparticle,yparticle,tn,savedir)
end
%
