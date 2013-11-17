

tn=10000

    
    
tn=tn+1;
    tic
    % findig nodes which are in the grain boundaries and solve
    % differential equation only for that points.
    [yii,xjj]=find(...
        imerode((phi>0.999),se)==0);
    % space discretization loop
    for ii=1:length(xjj)
        i=yii(ii);j=xjj(ii);
        % finding etas that are not zero
        etaind=find((eta(i,j,:)>0.001));
        del2=zeros(1,1,p);
        sumterm=del2;
        detadtM=del2;
        detadt=del2;

        del2(1,1,etaind)=1/delx^2*(0.5*(eta(indg(i+1,nboxsize),j,etaind)-2*eta(i,j,etaind)+eta(indg(i-1,nboxsize),j,etaind))...
            +0.25*(eta(indg(i+2,nboxsize),j,etaind)-2*eta(i,j,etaind)+eta(indg(i-2,nboxsize),j,etaind)))...
            +1/delx^2*(0.5*(eta(i,indg(j+1,mboxsize),etaind)-2*eta(i,j,etaind)+eta(i,indg(j-1,mboxsize),etaind))...
            +0.25*(eta(i,indg(j+2,mboxsize),etaind)-2*eta(i,j,etaind)+eta(i,indg(j-2,mboxsize),etaind)));
        sumterm=eta(i,j,etaind)*sum(eta(i,j,etaind).^2)-eta(i,j,etaind).^3;
        detadtM=(-alpha*eta(i,j,etaind)+beta*eta(i,j,etaind).^3-kappa*del2(1,1,etaind)+...
            2*epsilon*eta(i,j,etaind)*ppf(i,j));
        detadt=-L*(detadtM+2*gamma*(sumterm));
        eta2(i,j,etaind)=eta(i,j,etaind)+delt(tn)*detadt;

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

    % display speed of this step calulation. Using
    toc
    