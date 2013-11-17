%% 3D slice sparce reading
figure
clear
tn=	4000;
mboxsize=200;
nboxsize=200;
lboxsize=200;
delx=2;
x=[0:nboxsize-1]*delx;
y=[0:mboxsize-1]*delx;
z=[0:lboxsize-1]*delx;
[x,y,z]=meshgrid(x,y,z);
savedir='/media/disk/sim_res/3D_Nparticle/1/'
phidata=importdata([savedir 'result_ ' num2str(tn) ' .txt']);
phi=ones(mboxsize,nboxsize,lboxsize);
for n=1:size(phidata,1)
    phi(phidata(n,1)+1,phidata(n,2)+1,phidata(n,3)+1)=phidata(n,4);
end
slice(x,y,z,phi,[0 mboxsize*delx],[0 nboxsize*delx],[0 lboxsize*delx]);
axis([0 mboxsize*delx+1 0 nboxsize*delx+1 0 lboxsize*delx+1])
colormap gray
shading interp
axis equal
box on
colorbar
h=gcf;
h1=get(h,'children');
set(h1(1),'clim',[0.5 1])

%% DelG vs velocity
L=1;
m=1;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
delG=[.3 .2 .1 0.05];
dVdt=[1.874e4 4472 1834 603];
veloc=dVdt/(200*2)^2;
figure
plot(delG/intenergy,veloc/intenergy/mobility,'s-')
xlabel('Driving Pressure/\sigma_{gb}')
ylabel('Interface Velocity/M/\sigma_{gb}')
grid on
hold on
plot([0 0.3],[0 0.3],'g')


