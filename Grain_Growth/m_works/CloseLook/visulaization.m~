% return
figure
graymap=255/(max(max(phi-min(min(phi)))))*(phi-min(min(phi)));
%         subplot(2,1,1)
imshow(uint8(graymap));
%         subplot(2,1,2)
%         plot(reshape(eta(:,4,:),gridn,p))
title(strcat('Time= ', num2str(tn)))
figure(gcf);


%Movie
clear
figure
ni=0
for tn=10:5:2000
    filename=strcat('/home/cenna/Documents/Projects/Phase Transformation/Phase Field Model/Grain Growth/works/BMdata/bm',num2str(tn),'.mat');
    load(filename)
    ni=ni+1;
    graymap=255/(max(max(phi-min(min(phi)))))*(phi-min(min(phi)));
    %         subplot(2,1,1)
    imshow(uint8(graymap));
    %         subplot(2,1,2)
    %         plot(reshape(eta(:,4,:),gridn,p))
    title(strcat('Time= ', num2str(tn)))
%     M(ni)=getframe

end
%     imshow(uint8(graymap));
movie(M,20)
movie2avi(M,'uniform.avi')

%% Animation 3D
clear
dirstring=[pwd '/test02/'];
load(strcat(dirstring,'settings','.mat'))
figure
ni=0;
i=1;
Mtn=[1:3000];
% Mtn=[10000:2:10100 10100:100:20000 20200:200:40000];
% Mtn=[1480]
for tn=620:10:850
%     figure
    filename=strcat(dirstring,'/',num2str(tn),'.mat');
    load(filename)
    
    % drawing phi field
%     phi=sum(eta(:,:,:,1:p).^2,4);
%     drawgrains3D(phi,xparticle,yparticle,zparticle,tn)
    % Map energy field
    [ME,E]=calculateE_3D(eta,ppf,delx,settings);
    drawgrains3D(ME,xparticle,yparticle,zparticle,tn)
    MME(tn)=E;
    pause(0.03)
    % M(i)=getframe;
%         filename=strcat(pwd,'/',num2str(tn),'.png');
%         print('-f1','-dpng','-r300',filename)
end


% Energy contour plots


% energy  change plots
tn=[1:5:3000]*delt
figure
plot(tn,Epart)

%SUBPLOTS
clear all
costumstring='onefield3D_e2'
load(strcat(pwd,'/',costumstring,'/','setings','.mat'))
figure
ni=1
for tn= [50 200 400 600 800 850]
    filename=strcat(pwd,'/',costumstring,'/',num2str(tn),'.mat');
    load(filename);
    %     phi=sum(eta(:,:,1:p).^2,3);
    %     phi=phi+ppf;
    %      drawgrains(phi,xparticle,yparticle,tn,eta,ppf)
    phi=sum(eta(:,:,:,1:p).^2,4);
    phi=phi+ppf;
    %  drawgrains3D(phi,xparticle,yparticle,zparticle,tn)
    subplot(3,2,ni)
    drawisosurf(eta(:,:,:,2),xparticle,yparticle,zparticle,tn)
    %     contourgrains(eta,xparticle,yparticle,tn)
    %    drawE(rot90(rot90(E)),xparticle,yparticle,tn,eta,ppf)
    pause(0.01)
    %     Epart(ni)=particle_energy(E,xparticle,yparticle,(2+4)*scale,mboxsize,nboxsize);
    ni=ni+1;
    view([61 12])
    hold off
end


% eta profiles

figure
hold on
xline=1
eta1=eta(:,xline,1);
eta2=eta(:,xline,2);
x=[1:mboxsize]*delx;
plot(x,eta1,'b')
hold on
plot(x,eta2,'r')
plot(x,ppf(:,end/2),'g')


%--------------------------------------------------------------------------



clear

Mcostumstring=['scale4D15_smalle'; 'scale4D15_0.5e  '; 'scale4D15_1e    ';...
    'scale4D15_2e    '; 'scale4D15_4e    ' ;'scale4D15_5e    '];
Mcostumstring=cellstr(Mcostumstring)
for cc=1:6

    load(strcat(pwd,'/',Mcostumstring{cc},'/','setings','.mat'))
    figure
    tn=200;
    filename=strcat(pwd,'/',Mcostumstring{cc},'/',num2str(tn),'.mat');
    load(filename);
    phi=sum(eta(:,:,1:p).^2,3);
    phi=phi+ppf;
    %     drawgrains(phi,xparticle,yparticle,tn,eta,ppf)
    contourgrains(eta,xparticle,yparticle,tn*delt1)
end

% -------------------------------------------------------------------------
%% Stationary Straight Grain Boundary Energy
% 1. One Phase field with particle in the middle
clear all
global nboxsize mboxsize
global delx
scale=5;
x=linspace(0.25,0.75,10);
epsilon=5;
ni=1
for ni=1:length(x)
    [eta,ppf,ME,E]=stationarymin_p1(x(ni));
    Meta(:,:,ni)=eta;
    MEM(:,:,ni)=ME;
    drawE(ME,x(ni),eta,ppf)
    MinterfaceE(ni)=E/nboxsize/delx;
end
for ni=1:length(x)
    ME=MEM(:,:,ni);
    eta=Meta(:,:,ni);
    drawE(ME,x(ni),eta,ppf)
end
figure
plot(x,MinterfaceE,x,MinterfaceE,'.')


%% ---- ENERGY WELL ---- 
%% LOOP for energy of the system with one particle and interface
% 2. A grain boundary with particle in the middle
clear all
global nboxsize mboxsize
global delx
scale=4;
epsilon=linspace(1,10,5)
for ieps=1:20
    x=linspace(0.2,0.8,20);
    for ni=1:length(x)
        [eta,ppf,ME,E]=stationarymin_p2(x(ni),epsilon(ieps));
        MEM(:,:,ni)=ME;
        MinterfaceE(ni)=E/(nboxsize*delx);
    end
    WellProfile(ieps,:)=MinterfaceE;
    % measuring well depth
    % interpolation of well
    Wellx=linspace(0,mboxsize*delx,1000);
    Welly=interp1(fix(mboxsize*x)*delx,WellProfile,Wellx,'cubic');
    WellDepth(ieps)=Minterface(1)-min(Welly);
    
end

figure
plot(fix(mboxsize*x)*delx,WellProfile)
% -------------------------------------------------------------------------
%% --- looping for the 2 independent variables of epsilon and particle size
%%%% SEE run_energy_wells.m
% -------------------------------------------------------------------------

% potting energy wells for particle radius ipr
figure
ipr=1
plot(fix(mboxsize*x)*delx,...
    reshape(WellProfile(:,ipr,:),size(WellProfile,1),size(WellProfile,3)))
title(['Potentia Well For Particle Radius=' num2str(Pr(ipr))])

%%%% Mesh: Energy Well Depth
pr=linspace(5,1,5);
epsilon=linspace(1e-18,10e-18,5);

[Pr,EPS]=meshgrid(pr,epsilon);

pri=linspace(max(pr),min(pr),40);
epsiloni=linspace(min(epsilon),max(epsilon),50);

[Pri,EPSi]=meshgrid(pri,epsiloni);
WellDepthi=interp2(Pr,EPS,WellDepth,Pri,EPSi,'cubic');

figure
mesh(Pri,EPSi,WellDepthi)
hold on
mesh(Pr,EPS,WellDepth)
xlabel('Particle Radius (nm)')
ylabel('\epsilon Parameter (J/m^3)')
zlabel('Energy Well Depth')

%% 2D PLOT
pr=linspace(5,1,5);
epsilon=linspace(1e-18,10e-18,5);
[Pr,EPS]=meshgrid(pr,epsilon);

pri=linspace(max(pr),min(pr),40);
epsiloni=linspace(min(epsilon),max(epsilon),5);
[Pri,EPSi]=meshgrid(pri,epsiloni);
WellDepthi=interp2(Pr,EPS,WellDepth,Pri,EPSi,'cubic');

figure
plot(pr,WellDepth,'o')
hold on
plot(pri,WellDepthi)

%%%% Mesh: Energy Well Width
figure
pr=linspace(5,1,10);
epsilon=linspace(1e-18,10e-18,5);
[Pr,EPS]=meshgrid(pr,epsilon);

pri=linspace(max(pr),min(pr),40);
epsiloni=linspace(min(epsilon),max(epsilon),50);
[Pri,EPSi]=meshgrid(pri,epsiloni);

WellWidthi=interp2(Pr,EPS,WellWidth,Pri,EPSi,'cubic');

mesh(Pr,EPS,WellWidth)
hold on
mesh(Pri,EPSi,WellWidthi)
contour3(Pri,EPSi,WellWidthi,5)
xlabel('Particle Radius (nm)')
ylabel('\epsilon Parameter (J/m^3)')
zlabel('Energy Well Width')

%% 2D PLOT
pr=linspace(5,1,10);
epsilon=linspace(1e-18,10e-18,5);
[Pr,EPS]=meshgrid(pr,epsilon);
pri=linspace(max(pr),min(pr),40);

epsiloni=linspace(min(epsilon),max(epsilon),5);
[Pri,EPSi]=meshgrid(pri,epsiloni);
WellWidthi=interp2(Pr,EPS,WellWidth,Pri,EPSi,'cubic');

figure
plot(pr,WellWidth,'o')
hold on
plot(pri,WellWidthi)
xlabel('Particle Radius (nm)')

ylabel('Energy Well Width')



%% changing BC value and tracking particle-matrix interface energy
MBCValue=linspace(0,1,10)
for i=1:10
    [eta,ppf,ME,E]=stationarymin_p2(0.1,MBCValue(i));
    PME(i)=E;
end


% Analytical value for the energy of one phase field:
L=[1];
alpha=[1];
beta=[1];
gamma=1;
kappa=[1];
minE=1/4*alpha*(2*beta-alpha)/beta^2

syms eta
f0=-alpha/2*eta^2+beta/4*eta^4+minE

sigma=int(1*sqrt(2*kappa*f0),'eta',-1/beta*(beta*alpha)^(1/2),+1/beta*(beta*alpha)^(1/2))
sigma=eval(sigma)

% equilibrium profile
etai=linspace(1,-1,1000);
for i=1:length(etai)
    xi(i)=int(sqrt(kappa/2/f0),0,etai(i));
end
figure
xi=eval(xi);
xi=sort(xi);
plot((xi),-etai)

x=int(sqrt(kappa/2/f0))
etai=linspace(-0.9999,0.9999,1000);
xi=subs(x,'eta',etai);
% interface thickness
disorder=(1-etai.^2);
plot(xi,disorder)
wi=trapz(xi,disorder)

% analytical interface profile and thickness
xi=int(sqrt(kappa/2/f0))
xx

% Gerneralize solution

syms eta alpha beta
E=-alpha/2*eta^2+beta/4*eta^4;
extre=solve(diff(E));
minval=subs(E,'eta',extre(2));
E=E+minval;

alpha=[1];
beta=[1];
E=eval(E)
sigma=int(eval(2*sqrt(kappa*E)),'eta',eval(extre(3)),eval(extre(2)))
eval(sigma)

% Analytical attempt to solve energy of GB
L=[1];
alpha=[1];
beta=[1];
gamma=1;
kappa=[1];
syms eta1 eta2
E=-alpha/2*eta1.^2+beta/4*eta1.^4-alpha/2*eta2.^2+beta/4*eta2.^4 ...
    +gamma*eta1.^2*eta2.^2;

dE=diff(subs(E,eta1,0),'eta2')
extre=solve(dE)
minE=subs(E,'eta1',extre(2));
minE=1/4*alpha*(2*beta-alpha)/beta^2;
E=E+minE;

%% local free energy density of homogeneous System
%% 1. One phase filed:
L=[1];
alpha=[1];
beta=[1];
gamma=1;
kappa=[2];
minE=1/4*alpha*(2*beta-alpha)/beta^2;
eta=linspace(-2,2,1000);
E=-alpha/2*eta.^2+beta/4*eta.^4+minE;
figure;
plot(eta,E)

%% 2. Two phase filed
clear

m=2;
gamma=1.5;

eta1=linspace(-0.1,1.1,80);
eta2=eta1;
[Eta1,Eta2]=meshgrid(eta1,eta2);
E=m*(-1/2*Eta1.^2+1/4*Eta1.^4-1/2*Eta2.^2+1/4*Eta2.^4+gamma*(Eta1.^2.*Eta2.^2+Eta1.^2.*Eta2.^2)+0.25);
% with additional driving pressure:
G1=0.0;G2=0;
DelG12=G2-G1; DelG21=G1-G2;
phi12=(-Eta1+Eta2+1)/2;  phi21=(-Eta2+Eta1+1)/2; 
E=m*(-1/2*Eta1.^2+1/4*Eta1.^4-1/2*Eta2.^2+1/4*Eta2.^4+gamma*(Eta1.^2.*Eta2.^2+Eta1.^2.*Eta2.^2)+0.25)...
    +3*(phi12.^2/2-phi12.^3/3)*DelG12+3*(phi21.^2/2-phi21.^3/3)*DelG21+(G1+G2)/2;

figure
mesh(Eta1,Eta2,E)
xlabel('\eta_1')
ylabel('\eta_2')
zlabel('Energy Density')
h=get(gcf,'children');
axis([-0.1 1.1 -0.1 1.1 0 1])
set(h,'clim',[0 0.5])

%% drawing energy profile from equilibrium profile 
m=2
kappa=2
% Eta1=eta(1,:);
% Eta2=eta(2,:);
x=linspace(0,35,1000);
Eta1=0.5*(1-tanh(m/2/kappa*(x-17.5)));
Eta2=0.5*(1+tanh(m/2/kappa*(x-17.5)));

ME=m*(-1/2*Eta1.^2+1/4*Eta1.^4-1/2*Eta2.^2+1/4*Eta2.^4+1.5*(Eta1.^2.*Eta2.^2)+0.25)+...
    kappa/2*(gradient(Eta1).^2+gradient(Eta2).^2);

figure
plot(x,ME)
trapz(x,ME)
intenergy=1/3*sqrt(2*m*kappa)
%% with additional driving pressure:
G1=0.0;G2=0.0;
DelG12=G2-G1; DelG21=G1-G2;
phi12=(-Eta1+Eta2+1)/2;  phi21=(-Eta2+Eta1+1)/2; 
ME=m*(-1/2*Eta1.^2+1/4*Eta1.^4-1/2*Eta2.^2+1/4*Eta2.^4+gamma*(Eta1.^2.*Eta2.^2+Eta1.^2.*Eta2.^2)+0.25)...
    +3*(phi12.^2/2-phi12.^3/3)*DelG12+3*(phi21.^2/2-phi21.^3/3)*DelG21+(G1+G2)/2;

hold on
plot3(Eta1,Eta2,ME,'-o')
% plot3(eta(:,fix(nboxsize/2),1),eta(:,fix(nboxsize/2),2),ME(:,fix(nboxsize/2)),'g')






%% particle and onephase field interaction:
clear all
global nboxsize mboxsize
global delx
scale=5;
epsilon=linspace(0,20,21)
for ieps=1:21
    ni=ieps;
    [eta,ME,E]=particle_onefield_1D(epsilon(ieps));
    Meta(:,:,ni)=eta;
    MEM(:,:,ni)=ME;
    [ppf]=zeros(mboxsize,nboxsize);
    ppf(1:fix(mboxsize/2),:)=1;
    drawE(ME,0,0,epsilon(ieps),eta,ppf)
    MinterfaceE(ni)=E/nboxsize/delx;

end

figure
plot(epsilon,MinterfaceE)

%% local energy density of particle one phase field
clear all
global nboxsize mboxsize
global delx
scale=4;

alpha=[1];
beta=[1];
epsilon=5;
kappa=[2];
eta1=linspace(-1.5,1.5,80);
ppf=eta1;
[Eta1,PPF]=meshgrid(eta1,ppf);
E=-alpha/2*Eta1.^2+beta/4*Eta1.^4+epsilon*Eta1.^2.*PPF.^2;
minE=1/4*alpha*(2*beta-alpha)/beta^2;
E=E+minE;
figure
mesh(Eta1,PPF,E)
shading interp
colorbar
xlabel('\eta_1')
ylabel('\epsilon')
zlabel('Energy Density')
h=get(gcf,'children');
set(h(2),'clim',[0 0.5])
axis([-1.5 1.5 -1.5 1.5 0 1])

% adding profile from equilibrium calculations to the local energy surface
[eta,ME,E]=particle_onefield(epsilon);
[ppf]=zeros(mboxsize,nboxsize);
ppf(1:fix(mboxsize/2),:)=1;
xpos=fix(nboxsize/2);
eta1=eta(:,xpos,1);
ppf1=ppf(:,xpos);
ei=-alpha/2*eta1.^2+beta/4*eta1.^4+epsilon*eta1.^2.*ppf1.^2+minE;
hold on
plot3(eta1,ppf1,ei)
eta1=[eta(1,xpos,1) eta(end,xpos,1)];
ppf1=[ppf(1,xpos) ppf(end,xpos)];
ei=-alpha/2*eta1.^2+beta/4*eta1.^4+epsilon*eta1.^2.*ppf1.^2+minE;
plot3(eta1,ppf1,ei,'o')
figure
drawE(ME,0,0,epsilon,eta,ppf)


%% Velocity fields maps
h=subplot(2,3,2);
drawvelocity(Mvelocx);
set(h,'clim',[-0.5 0.5])
title('X velocity component')

h=subplot(2,3,3);
drawvelocity(Mvelocy)
set(h,'clim',[-0.7 0.7])
title('Y velocity component')

h=subplot(2,3,5);
drawvelocity(Mveloc)
set(h,'clim',[-0.9 0])
title('Velocity Magnitude')
xlabel(strcat('Time= ', num2str(tn*delt)))

h=subplot(2,3,6);
drawvelocity(nablaetax)
set(h,'clim',[-0.3 0.3])
title('\nabla \eta x')

h=subplot(2,3,4);
drawvelocity(nablaetay)
set(h,'clim',[-0.3 0.3])
title('\nabla \eta y')

h=subplot(2,3,1);
drawvelocity(Mdetadt)
set(h,'clim',[-0.2 0])
title('d\eta / dt')



%% Particles Drag

% local free energy density of homogeneous System
%% 1. One phase filed:
L=[1];
alpha=[-1];
beta=[1];
gamma=1;
nu1=0.5;
nu2=5;
kappa=[2];
minE=1/4*alpha*(2*beta-alpha)/beta^2;
eta=linspace(-1.5,1.5,1000);
E=-alpha/2*eta.^2+beta/4*eta.^4+minE-nu1*exp(-nu2*eta.^2);
figure;
plot(eta,E)
grid

%% 2. Two phase filed
alpha=[1];
beta=[1];
gamma=2;
nu1=0.2;
nu2=10;
kappa=[1];
eta1=linspace(-1.5,1.5,80);
eta2=eta1;
[Eta1,Eta2]=meshgrid(eta1,eta2);
E=-alpha/2*Eta1.^2+beta/4*Eta1.^4-alpha/2*Eta2.^2+beta/4*Eta2.^4 ...
    +gamma*Eta1.^2.*Eta2.^2-nu1*exp(-nu2*Eta1.^2)-nu1*exp(-nu2*Eta2.^2);
minE=1/4*alpha*(2*beta-alpha)/beta^2;
E=E+minE;
figure
mesh(Eta1,Eta2,E)
xlabel('\eta_1')
ylabel('\eta_2')
zlabel('Energy Density')
h=get(gcf,'children');
axis([-1.5 1.5 -1.5 1.5 0 1])
set(h,'clim',[0 0.5])



%% calculating gradient^2 of one phase field:
for i=1:mboxsize
    for j=1:nboxsize
        Del2(i,j)=1/delx^2*(0.5*(eta(indg(i+1,mboxsize),j)-2*eta(i,j)+eta(indg(i-1,mboxsize),j))...
            +0.25*(eta(indg(i+2,mboxsize),j)-2*eta(i,j)+eta(indg(i-2,mboxsize),j)))...
            +1/delx^2*(0.5*(eta(i,indg(j+1,nboxsize))-2*eta(i,j)+eta(i,indg(j-1,nboxsize)))...
            +0.25*(eta(i,indg(j+2,nboxsize))-2*eta(i,j)+eta(i,indg(j-2,nboxsize))));
%         detadtM=(-alpha.*eta(i,j)+beta.*eta(i,j).^3-kappa.*del2+...
%             2*epsilon.*eta(i,j)*ppf(i,j).^2);
        DetadtM(i,j)=(-alpha.*eta(i,j)+beta.*eta(i,j).^3+...
            +nu1*nu2*eta(i,j).*exp(-nu2*eta(i,j).^2)+...
            2*epsilon.*eta(i,j)*ppf(i,j).^2)...
            -kappa.*Del2(i,j);
        Detadt(i,j)=-L.*(DetadtM(i,j));
    end
end


%% Prossessing snapshots
dirstring='/media/disk/sim_res/triple_particle01';
% dirstring='/Drive2/sim_res/Particles2/results/';
load(strcat(dirstring,'/','settings','.mat'))

tn=4500
figure
filename=strcat(dirstring,'/',num2str(tn),'.mat');
load(filename)
drawvelocity(ME)

%% ------------------------------------------------------------------------
%% Energy Shift in Chen and Fan functional
%% ------------------------------------------------------------------------

% calculating velocity of the interface

timevec=[0:1:timestepn-1]*delt
figure
plot(timevec,MetaVol)


%% potential well in dynamic calculations


% ------ Dynamic Interface -------

figure
width=nboxsize*delx;
eqpos=MetaVol/width;
eqpos=Mintpos;
plot(eqpos,MME)
hold on
%fitting two straight lines on energy curve
ln1ind=find(Mintpos>16 & Mintpos<18);
plot(eqpos(ln1ind),MME(ln1ind),'y')
pp1=polyfit(eqpos(ln1ind),MME(ln1ind),1);

ln2ind=find(Mintpos>45);
plot(eqpos(ln2ind),MME(ln2ind),'g')
pp2=polyfit(eqpos(ln2ind),MME(ln2ind),1);

% potential well re-construction

well=MME-polyval((1*pp1+0*pp2),eqpos);
% well=MME-polyval(pp1,MetaVol/width);

figure
plot(MetaVol/width,well)
xlabel('Equivalent Interface Position (nm)')
ylabel('Total System Energy (J)')


% position vs. time
eqpos=Mintpos;
pf=sum(sum(ppf))/mboxsize/nboxsize;
t=timevec(2:end);
figure
plot(t,eqpos,'r')
% interface speed
pp=polyfit(t,eqpos,1);
intvel=pp(1)

% Average velocity
vaveg=(eqpos(end)-eqpos(10))/(t(end)-t(10))

%finding pinning force

intvel0=1.3396;
Pf=DelG(2)*(1-intvel/intvel0)
VelocityFromModelParameters=3/4*L(1)*sqrt(8*kappa(1)/alpha(1))*DelG(2)
% fraction of particles
pf=sum(sum(ppf))/(xend-x)/nboxsize;
sigma=4*kappa(1)/3*sqrt(8*kappa(1)/alpha(1))
% sigma=1.9521
PzFromZenerTheory=3/4*pf*sigma/pr

