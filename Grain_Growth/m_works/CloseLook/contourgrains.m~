function contourgrains(eta,xparticle,yparticle,tn,ppf)
% For making white regions brighter
global mboxsize nboxsize
global delx delt
p=2;
phi=sum(eta(:,:,1:p).^2,3);
phi(phi>1)=nan;
graymap=255/(max(max(phi-min(min(phi)))))*(phi-min(min(phi)));
mag=1;

% imshow(imresize(uint8(graymap),mag));
imshow(uint8(graymap));
hold onlinear
plot(xparticle*mag,yparticle*mag,'r.')
title(strcat('Time= ', num2str(tn*delt)))


contour(((eta(:,:,2))),[0.5 0.5],'color','r')
contour(((eta(:,:,1))),[0.5 0.5],'color','b')

axis equal
axis([0 nboxsize 0 mboxsize])
hold off
