
figure
graymap=255/(max(max(phi-min(min(phi)))))*(phi-min(min(phi)));
%         subplot(2,1,1)
imshow(uint8(graymap));
%         subplot(2,1,2)
%         plot(reshape(eta(:,4,:),gridn,p))
title(strcat('Time= ', num2str(tn)))
figure(gcf);


%% Movie
clear
figure
ni=0
i=1;
for tn=50:5:2000
    filename=strcat('/media/disk/sim_res/set0/',num2str(tn),'.mat');
    load(filename)
    ni=ni+1;
    graymap=255/(max(max(phi-min(min(phi)))))*(phi-min(min(phi)));
    %         subplot(2,1,1)
    imshow(uint8(graymap));
    %         subplot(2,1,2)
    %         plot(reshape(eta(:,4,:),gridn,p))
    title(strcat('Time= ', num2str(tn)))

%% statistics
    [areaG,diamG,perimG]=grainstat(phi);
    A{i}=areaG;
    D{i}=diamG;
    P{i}=perimG;
    Ai=A{i};Di=D{i};Pi=P{i};
    % removing small dots
%     [removeindex]=find(Ai<exp(1));
%     Ai(removeindex)=[];
%     Di(removeindex)=[];
%     Pi(removeindex)=[];
    An(i)=mean(Ai);
    %   Dn(i)=sqrt(mean((Di).^2));
    Dn(i)=mean(Di);
    Pn(i)=mean(Pi);
    GBArea(i)=sum(Pi)/2;
    i=i+1
end
%     imshow(uint8(graymap));


%% Animation
clear
dirstring='/media/disk/sim_res/set4/';
% dirstring='/media/disk/sim_res/Particles2/results/';
% dirstring='/Drive


2/sim_res/Particles2/dissolve1460/'
% dirstring='/media/disk/sim_res/SS_Band/'
% dirstring='/media/disk/sim_res/Particles1/Dissolve800/'
load(strcat(dirstring,'settings','.mat'))
figure
ni=0;
i=1;
% Mtn=[ 105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000];
% Mtn=[10000:2:10100 10100:100:20000 20200:200:40000];
Mtn=[5:20:100 140:40:1000 1120:120:4000 4200:200:10000]
for tn=Mtn
    filename=strcat(dirstring,'/',num2str(tn),'.mat');
    load(filename)
    drawgrains(phi,ppf,xparticle,yparticle,tn)
%     drawgrains(phi,0,0,0,tn)
    pause(0.01)
    % M(i)=getframe;
     filename=strcat(dirstring,'/',num2str(tn),'.png');
    %     imwrite(imadjust(phi),filename,'png')
    print('-f1','-r200','-dpng',filename)
end


%% SUBPLOT figures
clear
dirstring='/Drive2/sim_res/Particles1/Dissolve800/'
dirstring='/Drive2/sim_res/non_uniform_par1/'
dirstring='/media/disk/sim_res/set5'
load(strcat(dirstring,'/','settings','.mat'))
figure
ni=0;
sub=100*[0.3 0.6 1 4 6 10].^2;
for subi=1:6
    subplot(2,3,subi)
    filename=strcat(dirstring,'/',num2str(sub(subi)),'.mat');
    load(filename)
    %     drawgrains(phi,ppf,xparticle,yparticle,sub(subi))
    %     figure
    drawgrains(phi,0,0,0,sub(subi))
end



%% Time History Grains Statistics
clear
dirstring='/home/cenna/Office/media/disk/sim_res/run_friction_solute_01/';
% dirstring='/media/disk/sim_res/set0/';
load(strcat(dirstring,'/','settings','.mat'))
% figure
ni=0;
i=1;
for tn=[105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000]
    filename=strcat(dirstring,'/',num2str(tn),'.mat');
    load(filename)
    [areaG,diamG,perimG]=grainstat(phi);
    A{i}=areaG;
    D{i}=diamG;
    P{i}=perimG;
    i=i+1
    %     M(i)=getframe;
end
save(strcat(dirstring,'/','Grains','.mat'))

clear
% dirstring='/media/disk/sim_res/Friction06'
dirstring='/home/cenna/Office/media/disk/sim_res/run_friction_solute_01/';
load(strcat(dirstring,'/','Grains','.mat'))
i=1;
tn=[105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000];
try
    for tni=tn
        Ai=A{i};Di=D{i};Pi=P{i};
        % removing small dots
        [removeindex]=find(Ai<exp(1));
        Ai(removeindex)=[];
        Di(removeindex)=[];
        Pi(removeindex)=[];
        An(i)=mean(Ai);
        %   Dn(i)=sqrt(mean((Di).^2));
        Dn(i)=mean(Di);
        Pn(i)=mean(Pi);
        GBArea(i)=sum(Pi)/2;
        i=i+1;
    end
catch
    disp('Warning: simulation was incomplete')
end
% making Timehistory
Timehistory(1)=0;
if length(delt)==1
    tn=savetimesteps;
    delt=ones(1,max(savetimesteps))*delt
end
for ti=1:length(delt)-1
    Timehistory(ti+1)=Timehistory(ti)+delt(ti);
end


hold on
plot(Timehistory(tn(1:size(A,2))),An)
title('Area of Grains')

figure
plot(Timehistory(tn(1:size(A,2))),Dn)
title('Diameter of Grains')


%% kinetics in two different region:
figure
plot(Timehistory(tn(1:size(A,2))),An,'.')
title('Area of Grains')
hold on
%find the time of the first change based in the timestep number
ind1=find(Timehistory(tn)==Timehistory(800));
ind1=ind1-1;

plot(Timehistory(tn(1:ind1)),An(1:ind1),'k')

ind2=ind1+2;
ind3=find(Timehistory(tn)==Timehistory(3000));

plot(Timehistory(tn(ind2:ind3)),An(ind2:ind3),'r')

pp1=polyfit(Timehistory(tn(1:ind1)),An(1:ind1),1)
pp2=polyfit(Timehistory(tn(ind2:ind3)),An(ind2:ind3),1)

%% -------------Distributions-------------- %%

clear
dirstring='/Drive2/sim_res/Particles1/Dissolve800/'
dirstring='/Drive2/sim_res/non_uniform_par1/'
dirstring='/media/disk/sim_res/set0/'

load(strcat(dirstring,'Grains','.mat'))
tn=[52:2:100 105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000];
% find the index of tn with particular step number
step=2000;
tni=find(tn==step)
Ai=A{tni};Di=D{tni};Pi=P{tni};

% load the image of grains
filename=strcat(dirstring,num2str(step),'.mat');
load(filename);

%% Re-analysis of image if you are not sure:
[areaG,diamG,perimG,L]=grainstat_simple(phi);
Ai=areaG;
Di=diamG;
Pi=perimG;

% removing small dots
[removeindex]=find(Ai<exp(0));
Ai(removeindex)=[];
Di(removeindex)=[];
Pi(removeindex)=[];

GrainsA=Ai;
% figure
% [n,A]=hist(log(GrainsA),80);
% n=n/sum(n);
[ndist,Adist]=simplehist(log(GrainsA));

figure
 bar(Adist,ndist)
title('Number density histogram')
xlabel('Ln(Grain Area)')
ylabel('Frequency')

figure
bar(Adist,ndist.*exp(Adist)/sum(ndist.*exp(Adist)))
title('Area density histogram')
xlabel('Ln(Grain Area)')
ylabel('Frequency')

%% Extract two size population density Based on number density

[A1,A2,Delu]=ExtractDist(Adist(:),ndist(:),GrainsA,'yes');


%% Time History Distribution Analysis

clear
dirstring='/Drive2/sim_res/Particles1/Dissolve800/'
dirstring='/Drive2/sim_res/set0/'
load(strcat(dirstring,'Grains','.mat'))
tn=[52:2:100 105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000];
tn=[200:50:2000];
% h = waitbar(0,'Analysing Distribution');
try
    for step=tn;
        tni=find(tn==step);
        Ai=A{tni};Di=D{tni};Pi=P{tni};

        % removing small dots
        [removeindex]=find(Ai<exp(0));
        Ai(removeindex)=[];
        Di(removeindex)=[];
        Pi(removeindex)=[];
        % average properties
        An(tni)=mean(Ai);
        %   Dn(i)=sqrt(mean((Di).^2));
        Dn(tni)=mean(Di);
        Pn(tni)=mean(Pi);
        GBArea(tni)=sum(Pi)/2;

        % distribution analysis
        GrainsA=Ai;
        [ndist,Adist]=simplehist(log(GrainsA));
        [A1,A2,Delu]=ExtractDist(Adist(:),ndist(:),GrainsA,'no');
        pause(0.01)
        MA1(tni)=A1;
        MA2(tni)=A2;
        MDelu(tni)=Delu;
        %         waitbar(step/10000)
    end
catch
    disp('Warning: simulation was incomplete')
end
% close(h)
% making Timehistory
Timehistory(1)=0;
if length(delt)~=1
    for ti=1:length(delt)-1
        Timehistory(ti+1)=Timehistory(ti)+delt(ti);
    end
else
    for ti=1:max(tn)
        Timehistory(ti+1)=Timehistory(ti)+delt;
    end
end



%% Some time steps only:

clear
dirstring='/Drive2/sim_res/Particles1/Dissolve800/'
dirstring='/Drive2/sim_res/set7/'
dirstring='/media/disk/sim_res/Particles3/results/'
% load(strcat(dirstring,'Grains','.mat'))
tn=[52:2:100 105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000];
ni=0
Mstep=[3000:20:3100];
for step=Mstep
    ni=ni+1;
    tni=find(tn==step);
    %     Ai=A{tni};Di=D{tni};Pi=P{tni};
    % Analysis of image
    filename=strcat(dirstring,num2str(step),'.mat');
    load(filename);
    [areaG,diamG,perimG,L]=grainstat(phi);
    Ai=areaG;
    Di=diamG;
    Pi=perimG;

    % removing small dots
    [removeindex]=find(Ai<exp(2));
    Ai(removeindex)=[];
    Di(removeindex)=[];
    Pi(removeindex)=[];
    % average properties
    An(ni)=mean(Ai);
    %   Dn(i)=sqrt(mean((Di).^2));
    Dn(ni)=mean(Di);
    Pn(ni)=mean(Pi);
    GBArea(ni)=sum(Pi)/2;

    % distribution analysis
    GrainsA=Ai;
    [ndist,Adist]=simplehist(log(GrainsA));
    figure
    [A1,A2,Delu]=ExtractDist(Adist(:),ndist(:),GrainsA,'yes');
    title(strcat('tn= ',num2str(step)));
    MA1(ni)=A1;
    MA2(ni)=A2;
    MDelu(ni)=Delu;
end
figure
plot(MDelu,MA2./MA1,MDelu,MA2./MA1,'o')
figure
plot(Mstep,MDelu)

h=gcf
hch=get(gcf,'children')
set(hch(2),'LineWidth',1)
xlabel('Ln(Grain Area)')
ylabel('Number Density Probability Distribution')
% print -r300 -dpng -f1 'Particles2_dist_1460.png'


%% ***** Spacial distribution of particles and triple junctions ****
% with the function findtripple one can find position of tripple points in
% the field and then compare with position of the particles.

%% for plot
clear
dirstring='/media/disk/sim_res/triple_02/';
load(strcat(dirstring,'settings','.mat'))
Mtn=[ 50:2:100 105:5:1000 1020:20:4000 4050:50:10000 10100:100:20000 20200:200:40000];
tn=Mtn(200)
particles=[xparticle;yparticle];
figure
filename=strcat(dirstring,'/',num2str(tn),'.mat');
load(filename)
% make particles dark
phi=phi.*imcomplement(ppf);
drawgrains(phi,ppf,xparticle,yparticle,tn)
%     drawgrains(phi,0,0,0,tn)
[triples,BW]=findtripple(phi);
hold on
plot(triples(1,:),triples(2,:),'.','color','y')
%    find fraction of triple junctions occupied by particles
[tripleind,particleind]=TJonP(triples,particles);
plot(particles(1,particleind),particles(2,particleind),'o','MarkerSize',8)
plot(triples(1,tripleind),triples(2,tripleind),'yo','MarkerSize',8)
%superposition of BW image
figure
imshow(BW)
hold on
plot(xparticle,yparticle,'r.','MarkerSize',4)
plot(triples(1,:),triples(2,:),'.','color','y')
plot(particles(1,particleind),particles(2,particleind),'o','MarkerSize',8)
plot(triples(1,tripleind),triples(2,tripleind),'yo','MarkerSize',8)


%% for time history analysis
clear
dirstring='/media/disk/sim_res/triple_02/';
load(strcat(dirstring,'settings','.mat'))
Mtn=[ 50:2:100 105:5:1000 1020:20:4000 4050:50:10000];% 10100:100:20000 20200:200:40000];
particles=[xparticle;yparticle];
tni=0
for tn=Mtn(30:2:end)
    tni=tni+1;
    filename=strcat(dirstring,'/',num2str(tn),'.mat');
    load(filename)
    % make particles dark
    phi=phi.*imcomplement(ppf);
    [triples]=findtripple(phi);
    %    find fraction of triple junctions occupied by particles
    [tripleind,particleind]=TJonP(triples,particles);
    pinnedtriplef(tni)=length(tripleind)/length(triples);
    pinnedparticlef(tni)=length(particleind)/length(particles);
    totaltriplenum(tni)=length(triples);
    tnseq(tni)=tn; % storing sequence of tn that image processing is preforemd
    tn
end
% making Time history
Timehistory(1)=0;
if length(delt)~=1 % if delt is not constant
    for ti=1:length(delt)-1
        Timehistory(ti+1)=Timehistory(ti)+delt(ti);
    end
else  % if delt is constant
    for ti=1:max(tn) 
        Timehistory(ti+1)=Timehistory(ti)+delt;
    end
end

figure
plot(Timehistory(tnseq),pinnedtriplef)
title('Fraction of Pinned Triple Junctions')
xlabel('Time')
ylabel('f')
grid on

figure
plot(Timehistory(tnseq),pinnedparticlef)
title('Fraction of Particles at Tripple Junctions')
xlabel('Time')
ylabel('f')
grid on

figure
plot(Timehistory(tnseq),totaltriplenum)
title('Total Number of Tripple Junctions')
xlabel('Time')
ylabel('f')
grid on

figure
plot(Timehistory(tnseq),totaltriplenum.*pinnedtriplef)
title('Total Number of Pinned Junctions')
xlabel('Time')
ylabel('f')
grid on

figure
plot(Timehistory(tnseq),length(particles)./totaltriplenum.*pinnedtriplef)
title('Number of Particles/Number of triple junctions')
xlabel('Time')
ylabel('n/N*')
grid on



