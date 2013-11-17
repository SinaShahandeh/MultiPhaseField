%% Drawing gray scale map of the structure based on value of the phi field
%% particles are plotted by red markers
function drawgrains(phi,ppf,xparticle,yparticle,tn)
% For making white regions brighter because at particles position value of
% the phi is not completely 0 therefore at particles phi becomes larger
% than 1 and it decreases the contrast
phi(ppf==1)=min(min(phi));
% mapping value of phi to the range of the gray level
graymap=255/(max(max(phi-min(min(phi)))))*(phi-(min(min(phi))));

%% For grawing with more magnification
mag=1;
imshow(imadjust(imresize(uint8(graymap),mag)));
% 
% imshow(uint8(graymap));
title(strcat('Time= ', num2str(tn*0.25)))
hold on

plot(mag*xparticle,mag*yparticle,'r.','MarkerSize',4)
hold off
pause(0.1)


%% if you like colored one each color indicates one phase field. Add eta to
%% the input of this mfile and at the main program 
% colorimage=zeros(size(ppf));
% for i=1:size(eta,3)
% colorimage=eta(:,:,i)*i+colorimage;
% end
% imshow((phi.*colorimage)/max(max(phi.*colorimage)))
% colormap jet
% hold on
% 
% plot(mag*xparticle,mag*yparticle,'r.','MarkerSize',4)
% hold off
% pause(0.01)
% 
% colormap jet

