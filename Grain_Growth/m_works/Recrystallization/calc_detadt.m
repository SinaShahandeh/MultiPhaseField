%% Calculates the d(eta)/dt for the given phase field parameters
%% so we have a field consist of the energy geradient term if the field is
%% negative it means there is some change happening. The zero magnetude
%% means we have stationary interface or nothing at that point.

function [Mdetadt]=calc_detadt(eta,ppf,phi,setings)

L=setings.L;
alpha=setings.alpha;
beta=setings.beta;
gamma=setings.gamma;
kappa=setings.kappa;
epsilon=setings.epsilon;

se=strel('square',3);

[yii,xjj]=find(...
    imerode((phi>0.999),se)==0);
Mdetadt=zeros(size(phi));
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
    Mdetadt(i,j)=sum(reshape(-L*(detadtM+2*gamma*(sumterm)),1,size(sumterm,3)));
end
se=strel('disk',6);
Mdetadt=imcomplement(imdilate(ppf,se)).*Mdetadt;
%%ploting
figure
h=clf;
x=0:delx:(nboxsize-1)*delx;
y=0:delx:(mboxsize-1)*delx;
[X,Y]=meshgrid(x,y);

surf(X,Y,(((Mdetadt))))
set(get(h,'Children'),'clim',[-2e-4 2e-4])
title('d \eta / dt')
shading interp
view([0 90])
colormap jet
axis equal 
axis off
colorbar

