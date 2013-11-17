clear all
savedir='onefield';
mkdir(savedir)
% figure;
clf
% phase field parameters
L=3*[1 1];
alpha=[1 1];
beta=[1 1];
gamma=1;
kappa=[1 1];
epsilon=5;
G=[0 0];
% geometry settings
p=2;
global nboxsize mboxsize
global delx
scale=5;
mboxsize=16*scale; % y axis
nboxsize=25*scale; % x axis
delx=2/scale;
grainD=15*scale;

endtime=50;
timestepn=10000;
delt=endtime/timestepn;

%savesettings
save(strcat(pwd,'/',savedir,'/','setings.mat'))

% *** Phase Field Procedure *** (so small and simple piece of code!)
eta=zeros(mboxsize,nboxsize,p);
% making initial structure
eta(:,:,1)=circlegrain(mboxsize,nboxsize,nboxsize/2,3*mboxsize/4,grainD,'dome');
eta(:,:,2)=imcomplement(eta(:,:,1));
% ppf is phase variable representing particles
%particledistro(nboxsize,mboxsize,particles_number,radius)
[ppf,xparticle,yparticle]=particledistro(nboxsize,mboxsize,1,2*scale);

eta2=zeros(mboxsize,nboxsize,p); %pre-assignment
%initialization
for tn=1:1
    for i=1:mboxsize
        for j=1:nboxsize
            del2=1/delx^2*(0.5*(eta(indg(i+1,mboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-1,mboxsize),j,:))...
                +0.25*(eta(indg(i+2,mboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-2,mboxsize),j,:)))...
                +1/delx^2*(0.5*(eta(i,indg(j+1,nboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-1,nboxsize),:))...
                +0.25*(eta(i,indg(j+2,nboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-2,nboxsize),:)));
            sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
            detadtM=(-alpha.*reshape(eta(i,j,:),1,p)+beta.*reshape(eta(i,j,:),1,p).^3-kappa.*reshape(del2,1,p)+...
                2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j))+G;
            detadt=-L.*(detadtM+2*gamma*reshape(sumterm,1,p));
            eta2(i,j,:)=eta(i,j,:)+reshape(delt*detadt,1,1,2);
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
end

for tn=1:timestepn
    for i=1:mboxsize
        for j=1:nboxsize
            del2=1/delx^2*(0.5*(eta(indg(i+1,mboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-1,mboxsize),j,:))...
                +0.25*(eta(indg(i+2,mboxsize),j,:)-2*eta(i,j,:)+eta(indg(i-2,mboxsize),j,:)))...
                +1/delx^2*(0.5*(eta(i,indg(j+1,nboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-1,nboxsize),:))...
                +0.25*(eta(i,indg(j+2,nboxsize),:)-2*eta(i,j,:)+eta(i,indg(j-2,nboxsize),:)));
            sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
            detadtM=(-alpha.*reshape(eta(i,j,:),1,p)+beta.*reshape(eta(i,j,:),1,p).^3-kappa.*reshape(del2,1,p)+...
                2*epsilon.*reshape(eta(i,j,:),1,p)*ppf(i,j))+G;
            detadt=-L.*(detadtM+2*gamma*reshape(sumterm,1,p));
            eta2(i,j,:)=eta(i,j,:)+reshape(delt*detadt,1,1,2);
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
%     drawgrains(phi,xparticle,yparticle,tn,eta,ppf)
%     contourgrains(eta,xparticle,yparticle,tn,ppf)
[cur1,cur2,detach]=analysecontourgrains(eta,xparticle,yparticle,tn,ppf);
%     [deletax deletay]=gradient(eta,delx,delx);
%     gradeta=sqrt(deletax.^2+deletay.^2);
%     E=-alpha/2*eta.^2+beta/4*eta.^4+kappa*gradeta.^2;
%     ME(tn)=sum(sum(E))/mboxsize/nboxsize;
%     drawE(rot90(rot90(E)),xparticle,yparticle,tn,eta,ppf)
    %     savegrains(eta,ppf,E,xparticle,yparticle,tn,savedir)
curi1(tn)=cur1;
curi2(tn)=cur2;
detachi(tn)=detach;


subplot(3,1,3)
plot(curi1,'b');
hold on
plot(curi2,'r');
hold off
   title(strcat('\kappa_{\eta_1}= ', num2str(cur1),' , ',...
        '\kappa_{\eta_2}= ', num2str(cur2)));

pause(0.005)

end
%
