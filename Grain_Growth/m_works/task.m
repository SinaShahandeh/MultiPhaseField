clear
clf
p=36;
gridn=800;
delx=2;
delt=0.25;
L=1;
alpha=1;
beta=1;
gamma=1;
kappa=2;
eta=zeros(gridn,gridn,p);
for nn=1:35000;
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
    tn
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
    phi=sum(eta(:,:,1:p).^2,3);
    eta=eta2;
    if mod(tn,5)==0
        filename=strcat('C:\My docs\Projects\Phase Transformation\Phase Field Model\Grain Growth\works\',num2str(tn),'.mat')
        save(filename,'phi')
    end
        
end
%



