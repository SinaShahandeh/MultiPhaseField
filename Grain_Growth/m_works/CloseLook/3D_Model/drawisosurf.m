function drawisosurf(eta,xp,yp,zp,tn)
% For making white regions brighter
global mboxsize nboxsize lboxsize
global delx

% 
% x=[0:nboxsize-1]*delx;
% y=[0:mboxsize-1]*delx;
% z=[0:lboxsize-1]*delx;
% [x,y,z]=meshgrid(x,y,z);
% slice(x,y,z,eta,xp*delx,(mboxsize-1)*delx,zp*delx)
% shading interp
% colormap jet
% title(strcat('Time= ', num2str(tn)))
% hold on
% plot3(xp*delx,yp*delx,zp*delx,'bo')
% hold off
% colorbar
% pause(0.01)

x=[0:nboxsize-1]*delx;
y=[0:mboxsize-1]*delx;
z=[0:lboxsize-1]*delx;
[x,y,z]=meshgrid(x,y,z);
[f,v]=isosurface(x,y,z,eta,0.5); 
p=patch('Faces',f,'Vertices',v);
grid on
% axis vis3d;s
axis([0 nboxsize*delx 0 mboxsize*delx 0 lboxsize*delx])
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
% view(3);
view([110 20])
axis equal
camlight 
lighting gouraud
axis([0 nboxsize*delx 0 mboxsize*delx 0 lboxsize*delx])

title(strcat('Time= ', num2str(tn)))
xlabel('X')
ylabel('Y')
zlabel('Z')
hold on
plot3(xp*delx,yp*delx,zp*delx,'ko')
hold off
pause(0.1)
