return
clear
clf
p=4;
gridn=50;
delx=2;
delt=0.25;
L=1;
alpha=1;
beta=1;
gamma=1;
kappa=2;
eta=zeros(gridn,gridn,p);
for nn=1:30000;
    ii=fix(gridn*rand(1,1))+1;jj=fix(gridn*rand(1,1))+1;
    eta(ii,jj,fix(p*rand(1,1))+1)=1;
end

% eta(:,5,1)=1;
% eta(:,:,2)=1;
% eta(:,5,2)=0;
% eta=rand(gridn,gridn,p)*0.001;%-0.001;
eta2=zeros(gridn,gridn,p); %pre assignment
t=0:delt:500;
for tn=1:size(t,2)
    for i=1:gridn
        for j=1:gridn
            del2=1/delx^2*(0.5*(eta(indg(i+1),j,:)-2*eta(i,j,:)+eta(indg(i-1),j,:))...
                +0.25*(eta(indg(i+2),j,:)-2*eta(i,j,:)+eta(indg(i-2),j,:)))...
                +1/delx^2*(0.5*(eta(i,indg(j+1),:)-2*eta(i,j,:)+eta(i,indg(j-1),:))...
                +0.25*(eta(i,indg(j+2),:)-2*eta(i,j,:)+eta(i,indg(j-2),:)));
            sumterm=eta(i,j,:)*sum(eta(i,j,:).^2)-eta(i,j,:).^3;
            detadtM=(-alpha*eta(i,j,:)+beta*eta(i,j,:).^3-kappa*del2);     
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
        phi(:,:,tn)=sum(eta(:,:,1:p).^2,3);
        graymap=255/(max(max(phi(:,:,tn)-min(min(phi(:,:,tn))))))*(phi(:,:,tn)-min(min(phi(:,:,tn))));
%         subplot(2,1,1)
        imshow(uint8(graymap));
%         subplot(2,1,2)
%         plot(reshape(eta(:,4,:),gridn,p))
        title(strcat('Time= ', num2str(tn)))    
                figure(gcf);
    eta=eta2;
    tn
end
%


%% making one big grain and small nuclie arround it

eta(:,:,1)=0;
Bigcirc=108;
Bigcirc2=100;
for pn=1:p
    eta(nboxsize/2-Bigcirc/2+1:nboxsize/2+Bigcirc/2,mboxsize/2-Bigcirc/2+1:mboxsize/2+Bigcirc/2,pn)=...
        eta(nboxsize/2-Bigcirc/2+1:nboxsize/2+Bigcirc/2,mboxsize/2-Bigcirc/2+1:mboxsize/2+Bigcirc/2,pn)...
        .*imcomplement(imcircle(Bigcirc));
end
eta(nboxsize/2-Bigcirc2/2+1:nboxsize/2+Bigcirc2/2,mboxsize/2-Bigcirc2/2+1:mboxsize/2+Bigcirc2/2,1)=imcircle(Bigcirc2);

for repnum=1:120
    xt=[];yt=[];
    for chunki=1:40
        xchunk=2*rand(1,1)-0.5;
        ychunk=2*rand(1,1)-0.5;
        x=0.05*randn(1,300)+xchunk;
        y=0.05*randn(1,300)+ychunk;
        xt=[xt x];
        yt=[yt y];
    end
    [vxt,vyt]=voronoi(xt,yt);
    [vt,ct]=voronoin([xt' yt']);
    inti=find(xt<1 & xt>0 & yt<1 & yt>0);
    cint=ct(inti);
    for i=1:length(cint)
        Aint(i)=polyarea(vt(cint{i},1),vt(cint{i},2));
    end
    repA=[repA Aint];
end



%%---------------------------------------------------------------------
%% Recrystallization works
Mparam=linspace(0,100,10);
for ni=6:length(Mparam)
    savedir=['c:\PF_Data\test_guillaume\band_parN\' num2str(ni) '\']
    Mfraction{ni}=run_particles(Mparam(ni),savedir);
end

    
    