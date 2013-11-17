function drawvelocity(MV)
% For making white regions brighter
global mboxsize nboxsize
global delx


x=linspace(0,nboxsize*delx,nboxsize);
y=linspace(0,mboxsize*delx,mboxsize);
[X,Y]=meshgrid(x,y);

% h=surf(X,Y,rot90(rot90(MV)));
 h=surf(X,Y,((MV)));
axis([0 nboxsize*delx 0 mboxsize*delx -10 10])

% axis equal
view([0 90])
colormap jet
shading interp
% set(get(h,'parent'),'clim',[0 1.5]*1e-19)


colorbar
 axis equal
