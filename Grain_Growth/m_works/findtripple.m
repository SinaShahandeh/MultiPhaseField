% Finding tripple points. The tripples will be two row matrix first row for
% the x position of the tripple jucntion and 2nd row for the y position of
% the tripple junction. The unit of the result is in pixels of the finite
% difference grid.
function [tripples,BW]=findtripple(phi)

% imshow(phi)
phi(phi>1)=nan;
graymap=uint8(255/(max(max(phi-min(min(phi)))))*(phi-(min(min(phi)))));
%         subplot(2,1,1)
%         imshow(imadjust(imresize(uint8(graymap),mag),[10 255]./255,[]));
%         subplot(2,1,2)
%         plot(reshape(eta(:,4,:),gridn,p))
graymap=imadjust(graymap,[120 255]./255,[]);
grayind=graythresh(graymap);
BW=im2bw(graymap,grayind);
BW=imfill(BW,'holes');
BW=bwmorph(BW,'thicken',inf);
BW=imcomplement(bwmorph(imcomplement(BW),'thin',inf));
BW=imcomplement(bwmorph(imcomplement(BW),'spur',inf));

L = bwlabel(BW,4);
s  = regionprops(L, 'Area','Centroid');
% removing small objects
C=[];A=[];
for n=1:length(s)
    A(n)=s(n).Area;
    C(n,:)=s(n).Centroid;
end
[smallind]=find(A<16); % 15 is size of biggest small object
BWsmall=bwselect(BW,C(smallind,1),C(smallind,2),4);
BW=BW.*imcomplement(BWsmall);
BW=imcomplement(bwmorph(imcomplement(BW),'thin',inf));
BW=imcomplement(bwmorph(imcomplement(BW),'spur',inf));

% BW=imclearborder(BW,4);
%
%
% figure
% imshow(BW)
% hold on
%% detecting the boundaries, it turns out not a good way
% [B,L,N,A] = bwboundaries(BW,8);
%
% hold on
% for k = 1:length(B)
%     boundary = B{k};
%     plot(boundary(:,2), boundary(:,1), 'LineWidth', 1)
% end
% % making a matrix containing all the boundary points
% bgpoints=[];
% for k=1:length(B)
%     boundary=B{k};
%     bgpoints=[bgpoints; boundary];
% end
%

%% much simpler way
% fining all pixels at boundaries then making a box of 3x3 to check if
% there is 4 black pixels in it then its triple junction
[nboxsize,mboxsize]=size(BW);
% finding points in the bounudary
[boundariesi,boundariesj]=find(BW==0);
%
tripples=[];
for bi=1:length(boundariesj)
    i=boundariesi(bi); %row
    j=boundariesj(bi); %column
    boxm=BW(indg(i-1,nboxsize):indg(i+1,nboxsize),indg(j-1,mboxsize):indg(j+1,mboxsize));
    blackpix=find(boxm==0);
    if length(blackpix)==4
        addpoint=[j;i];
        tripples=[tripples addpoint];
    end
end
% some of the tripple points are repetetive in some locations for removing
% those we average over close points
% making distance matrix
Mdistance=zeros(length(tripples),length(tripples))+100;
for i=1:length(tripples)
    for j=1:i
        dx=tripples(2,i)-tripples(2,j);
        dy=tripples(1,i)-tripples(1,j);
        Mdistance(i,j)=sqrt(dx^2+dy^2);
    end
end
Mdistance=Mdistance+eye(length(tripples),length(tripples))*100;
% find close points
[iclose,jclose]=find(Mdistance<2);
% the final point is average of the two close points which is close enough
% to real triple point
tripples(:,jclose)=(tripples(:,iclose)+tripples(:,jclose))/2;
tripples(:,iclose)=[];

% plot(tripples(1,:),tripples(2,:),'*','color','y')


return
%% code check
figure
imshow(BW)
hold on
plot(tripples(1,:),tripples(2,:),'*','color','y')

