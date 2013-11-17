function drawgrains3D(phi,xp,yp,zp,tn)
% For making white regions brighter
global mboxsize nboxsize lboxsize
global delx
% phi(phi>1)=nan;

graymap=255/(max(max(max(phi-min(min(min(phi)))))))*(phi-min(min(min(phi))));
%         subplot(2,1,1)

x=[0:nboxsize-1]*delx;
y=[0:mboxsize-1]*delx;
z=[0:lboxsize-1]*delx;
[x,y,z]=meshgrid(x,y,z);
slice(x,y,z,graymap,xp*delx,[],zp*delx)
shading interp
colormap jet
title(strcat('Time= ', num2str(tn)))
% hold on
% plot3(xp*delx,yp*delx,zp*delx,'bo')
 view([90 0])
hold off
colorbar
pause(0.01)