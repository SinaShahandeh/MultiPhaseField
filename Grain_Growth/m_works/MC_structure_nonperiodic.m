
%% Read Image
adress='image.png'
phi=imread(adress);

figure
imshow(imadjust(phi))
% phi(phi>1)=nan;
graymap=uint8(255/(max(max(phi-min(min(phi)))))*(phi-(min(min(phi)))));

% detecting the grains and labeling the objects
graymap=imadjust(graymap);
grayind=graythresh(graymap);
BW=im2bw(graymap,grayind);
BW=imfill(BW,'holes');
BW=bwmorph(BW,'thicken',inf);
BW=imcomplement(bwmorph(imcomplement(BW),'thin',inf));
BW=imcomplement(bwmorph(imcomplement(BW),'spur',inf));
% BW=imclearborder(BW,4);
L = bwlabel(BW,4);

% repeat this loop untill there is empty a and b at the end
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


