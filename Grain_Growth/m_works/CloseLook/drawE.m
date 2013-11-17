function drawE(ME,tn,eta,ppf)
% For making white regions brighter
global mboxsize nboxsize
global delx
subplot(2,1,1)
x=0:delx:(nboxsize-1)*delx;
y=0:delx:(mboxsize-1)*delx;
[X,Y]=meshgrid(x,y);
surf(X,Y,ME)
shading interp
view([0 90])
axis off

colormap jet
axis equal 
colorbar
pause(0.003)
% imshow(uint8(graymap));

title(strcat('Time= ', num2str(tn)))
subplot(2,1,2)
xpos=fix(nboxsize/2);
eta1=eta(:,xpos,1);

x=[1:mboxsize]*delx;
plot(x,eta1,'r')
hold on
plot(x,ppf(:,xpos),'g')
hold off
pause(0.003)
