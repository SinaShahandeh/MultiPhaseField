function [areaG,diamG,perimG,L]=grainstat(phi)
% For making white regions brighter
[n,m]=size(phi);
phi(phi>1)=1;
graymap=uint8(255/(max(max(phi-min(min(phi)))))*(phi-(min(min(phi)))));
%         subplot(2,1,1) 
%         imshow(imadjust(imresize(uint8(graymap),mag),[10 255]./255,[]));
%         subplot(2,1,2)
%         plot(reshape(eta(:,4,:),gridn,p))
graymap=imadjust(graymap,[120 255]./255,[]);
grayind=graythresh(graymap);
BW=im2bw(graymap,grayind);
% removing particles inside the grains
BW=imfill(BW,'holes');

% finding the grains at border
BWborder=imclearborder(BW,4);
se=strel('square',2);
BWborderclose=imclose(BWborder,se); 
BWborderg=imcomplement(BWborder).*BW;

triBW=triu(BWborderg)+BWborder;
% putting Boundary grains beside the main box
% in periodic boundary condition
L = bwlabel(BWborderg,4);
s  = regionprops(L, 'BoundingBox');
s=struct2cell(s);
s=cell2mat(s);
% finding width and height of the objects at the border
sx=s(3:4:end);
sy=s(4:4:end);
dividx=max(sx)+2;
dividy=max(sy)+2;
BWbig=[BWborderg(n-dividy:end,:) BWborderg(n-dividy:end,1:dividx);...
    triBW BWborderg(:,1:dividx)];
BWbig=imclearborder(BWbig,4);
% making boundaries thin to one pixel width
BWbig=imcomplement(bwmorph(imcomplement(BWbig),'thin',3));
BWbig=imcomplement(bwmorph(imcomplement(BWbig),'spur',inf));
BWbig=imclearborder(BWbig,4);
% measuring objects properties
L = bwlabel(BWbig,4);
s  = regionprops(L, 'Area','EquivDiameter','MajorAxisLength',...
    'MinorAxisLength','Perimeter');
areaG = cat(1, s.Area);
% Diam=(cat(1, s.MajorAxisLength)+cat(1, s.MinorAxisLength))/2;
diamG=cat(1,s.EquivDiameter);
perimG=cat(1,s.Perimeter);
% 
% figure
% hold on
% plot(mag*xparticle,mag*yparticle,'r.','MarkerSize',4)
% hold off
% imshow(BW)
% pause(1)
% 
