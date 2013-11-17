return

%% Changing epsilon and traking particle-matrix interface energy

Mepsilon=[linspace(1,50,20) linspace(55,100,10)]*1e-18;
n=1
for epsilon=Mepsilon
    MSigma(n)=realtime_1D_particle_interface_func(epsilon)
    n=n+1
end

figure;
plot(Mepsilon,MSigma)
grid on
xlabel('\epsilon parameter')
ylabel('Particle-Matrix Interface Energy')

===========================================================================
%% +++++ Applying Driving Force on Interface ++++++
%% calculating velocity of 1D grain as a function of excess energy of the grain
clear all
Mparam=linspace(0,1.5e-19,10)
for nn=1:length(Mparam)
    MMetavol(nn,:)=interface_p2_1D(Mparam(nn))';
end

% times should be set by hand
timestepn=300;
delt=0.01;

% plotting volume of grain 1 vs. time
figure
timevec=[0:1:timestepn-1]*delt
plot(timevec,MMetavol)

% calculating velocity of the interface at different driving pressure
for nn=1:length(Mparam)
    ppvel=polyfit(timevec(200:300),MMetavol(nn,200:300),1);
    vel(nn)=ppvel(1);
end

figure
plot(Mparam,vel,'ro')
hold on

% plotting profiles and energy
figure
subplot(2,1,1)
hold on
eta1=eta(:,1);
eta2=eta(:,2);
x=[1:mboxsize]*delx-(MetaVol(end)-MetaVol(1));
plot(x,eta1,'b')
hold on
plot(x,eta2,'r')

subplot(2,1,2)
x=[1:mboxsize]*delx;
plot(x,ME(:,1),'r')

===========================================================================

%% potential well in dynamic calculations
timestepn=300;
delt=0.1;
% plotting energy of system vs. time
figure
timevec=[0:1:timestepn-1]*delt
plot(MME)

===========================================================================
%% === Friction of fine particles on the interface motion ===
%% changing Pz and track speed of the interface
clear
figure
savedir='/media/disk/sim_res/DomeDrivingForce_particle_3/';
savedir= [pwd '/']
Mparam=linspace(20,20,1);
for filenum=1:length(Mparam)
    param=Mparam(filenum)
    [InterfaceVel]=realtime_particle_friction(param,filenum,savedir,Pz);
    MInterfaceVel(filenum)=InterfaceVel;
%     Mpp(filenum)=pp
end
figure
plot(Mparam,MInterfaceVel)
grid

%interface mobiolity:
InterfaceCurvature=1./Mparam;
plot(1./Mparam,-MInterfaceVel)
sigma=4*kappa(1)/3*sqrt(8*kappa(1)/alpha(1))
pp=polyfit(1./Mparam,-MInterfaceVel,1)
Mobility=pp(1)/sigma
MobilityFromModelParameters=3/4*L(1)*sqrt(8*kappa(1)/alpha(1))
DrivingPressure=-MInterfaceVel./MobilityFromModelParameters;
constant=polyfit(InterfaceCurvature,DrivingPressure,1)

%% ploting Volume of dome vs time
% loading data from saved Eta Volume data
dirread='/media/disk/sim_res/DomeDrivingForce/'
figure
for filenum=1:20
    load([dirread num2str(filenum) '.mat'])
    hold on
    plot(timevec(2:end),MetaVol(1,:)/(mboxsize*nboxsize*delx^2))
    MInterfaceVel(filenum)=InterfaceVel;
end
Pz=[linspace(0,0.05,20)];
figure
plot(Pz,MInterfaceVel)
grid

===========================================================================
%% Changing Gamma and track energy of triple junction vs. gb

clear
MGamma=[linspace(0.5,5,40)];
for filenum=1:40
    Gamma=MGamma(filenum);
    [ETriple,EGB,E]=run_triple_one(0,Gamma,filenum)
    METriple(filenum)=ETriple;
    MEGB(filenum)=EGB;
    MEtot(filenum)=E;
end
save

figure
plot(MGamma,METriple./MEGB)

% reading energy data from saved simulation files
clear
savedir='/media/disk/sim_res/triple_variables01/';
MGamma=[linspace(0.5,5,40)];
for filenum=1:40
    load([savedir num2str(filenum) '.mat'])
    ind=findtripple(phi);
    ETriple=ME(fix(ind(2))+1,fix(ind(1))+1);
    EGB=max(ME(end,:));
    
    METriple(filenum)=ETriple;
    MEGB(filenum)=EGB;
    MEtot(filenum)=E;
end

figure
plot(MGamma,METriple./MEGB)


%% triple point with a particle on top of it


clear
MGamma=[linspace(0.5,5,40)];
for filenum=1:40
    Gamma=MGamma(filenum);
    [ETriple,EGB,E,Ep]=run_triple_one_particle(Gamma,filenum);
    METriple(filenum)=ETriple;
    MEGB(filenum)=EGB; 
    MEtot_N(filenum)=E;
    MEtot_P(filenum)=Ep;
end
save('triple_variable02.mat')

clear
load('triple_variable02.mat')

figure
plot(MGamma,METriple./MEGB)

figure
plot(MGamma,(MEtot_N-MEtot_P))


%% New functional for triple points
% Changing Gamma and track energy of triple junction vs. gb
clear
MGamma=[linspace(0.5,5,20)];
for filenum=1:40
    Gamma=MGamma(filenum);
    [ETriple,EGB,E]=run_triple_functional(Gamma,-1,filenum)
    METriple(filenum)=ETriple;
    MEGB(filenum)=EGB;
    MEtot(filenum)=E;
end

%retrive data from saved files
clear
MGamma=[linspace(0.5,5,20)];
for filenum=1:20
    Gamma=MGamma(filenum);
    load(['/media/disk/sim_res/triple_variables_functional_01/' num2str(filenum) '.mat'])
    ind=findtripple(phi);
    ETriple=ME(fix(ind(2))+1,fix(ind(1))+1);
    EGB=max(ME(end,:));
    ETriple/EGB
    
    METriple(filenum)=ETriple;
    MEGB(filenum)=EGB;
    MEtot(filenum)=E;
end

figure
plot(MGamma,METriple./MEGB)


%% Measuring Angle Manually

[x,y]=ginput(3);

v1=[x(1)-x(2) y(1)-y(2)];
v2=[x]






