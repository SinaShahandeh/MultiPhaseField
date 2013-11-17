
%% Read Image
adress='image.png'
phi=imread(adress);

figure
imshow(imadjust(phi))

% detecting the grains and labeling the objects
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

%% *manual* remove repetetive grains on the top corner of L

figure
imshow(L)

BW2=imfill(imcomplement(BWbig));
BWbig=imcomplement(BW2);
L = bwlabel(BWbig,4);

%% make labeled image (im) periodic based on the reconstructed image L
[n,m]=size(phi);
[i,j]=size(L);
im=L(i-n+1:end,1:m);
% change this 6 and 4 number based on the width of grain boundary (may cause
% problem and unaccuracy for some structures :(  )
im(:,1:j-m+1-6)=L(i-n+1:end,m+6:end)+im(:,1:j-m+1-6);
im(2*n-i+1+4:n,:)=L(1:i-n-4,1:m)+im(2*n-i+1+4:n,:);
im(2*n-i+1:n,1:j-m+1)=L(1:i-n,m:end)+im(2*n-i+1:n,1:j-m+1);

L=im;



%% *manual* repeat this loop untill there is empty a and b at the end
for j=1:size(phi,2)-1
    for i=2:size(phi,1)
        if L(i,j)==0
            L(i,j)=L(i,j+1);
        end
    end
end

for j=1:size(phi,2)
    for i=2:size(phi,1)
        if L(i,j)==0
            L(i,j)=L(i-1,j);
        end
    end
end

[a,b]=find(L==0);
% only upper row remaining
for j=1:size(phi,2)
    for i=1:size(phi,1)-1
        if L(i,j)==0
            L(i,j)=L(i+1,j);
        end
    end
end
[a,b]=find(L==0)


%% 

% converting labeled image to RGB color
RGBL=label2rgb(L,'jet','k','shuffle');
figure
imshow(RGBL)
RGBLd=double(RGBL);
% converting to orientation variable
ori=255*255*RGBLd(:,:,1)+255*RGBLd(:,:,2)+RGBLd(:,:,3);



% writer the structure into text file
fid=fopen('Structure_big.txt','wt')
for i=1:size(ori,2)
    for j=1:size(ori,1)
        fprintf(fid,'% 10.0f',ori(i,j));
    end
    fprintf(fid,'\n');
end

fclose(fid)


