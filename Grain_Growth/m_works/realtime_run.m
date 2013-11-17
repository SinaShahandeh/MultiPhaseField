clear
figure;
clf
p=2;
gridn=100;
delx=2;
delt=0.25;
L=1;
alpha=1;
beta=1;
gamma=1;
kappa=2;
eta=zeros(gridn,gridn,p);
% putting initial nucleas 1:total number of nuclea (so nucleation density
% can be set)
for nn=1:1;
%     ii=fix(gridn*rand(1,1))+1;jj=fix(gridn*rand(1,1))+1;
    ii=gridn/2;
    jj=gridn/2;
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
            del2=1/delx^2*(0.5*(eta(indg(i+1,gridn),j,:)-2*eta(i,j,:)+eta(indg(i-1,gridn),j,:))...
                +0.25*(eta(indg(i+2,gridn),j,:)-2*eta(i,j,:)+eta(indg(i-2,gridn),j,:)))...
                +1/delx^2*(0.5*(eta(i,indg(j+1,gridn),:)-2*eta(i,j,:)+eta(i,indg(j-1,gridn),:))...
                +0.25*(eta(i,indg(j+2,gridn),:)-2*eta(i,j,:)+eta(i,indg(j-2,gridn),:)));
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
    
    phi=sum(eta(:,:,1:p).^2,3);
    eta=eta2;


        graymap=255/(max(max(phi-min(min(phi)))))*(phi-min(min(phi)));
        %         subplot(2,1,1)
        imshow(imresize(uint8(graymap),3));
        %         subplot(2,1,2)
        %         plot(reshape(eta(:,4,:),gridn,p))
        title(strcat('Time= ', num2str(tn)))
%         M=getframe
pause(0.01)
end
%
