function drawg(phi,xparticle,yparticle,tn,eta,ppf)
% For making white regions brighter
global mboxsize nboxsize
global delx
phi(phi>1)=nan;

graymap=255/(max(max(phi-min(min(phi)))))*(phi-min(min(phi)));
%         subplot(2,1,1)

mag=1;
subplot(2,1,1)
imshow(imresize(uint8(graymap),mag));
% imshow(uint8(graymap));

title(strcat('Time= ', num2str(tn)))
hold on
plot(xparticle*mag,yparticle*mag,'r.')
contour(((eta(:,:,1))),[0 0],'color','y')
hold off
subplot(2,1,2)
xpos=fix(nboxsize/2);
eta1=eta(:,xpos,1);

x=[1:mboxsize]*delx;
plot(x,eta1,'r')
hold on
plot(x,ppf(:,xpos),'g')
hold off
pause(0.003)
