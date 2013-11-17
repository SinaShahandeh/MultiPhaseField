
%% Size Distribution  %%
% grain structure statistics at time step tn. Read the results from savedir directory
clear
savedir='Uniform'
load(strcat(pwd,'/',savedir,'/','setings','.mat'))
% figure
tn=2000
filename=strcat(pwd,'/',savedir,'/',num2str(tn),'.mat');
load(filename)
phi(logical(ppf))=0;
[areaG,diamG,perimG,L]=grainstat(phi);

[N,a]=hist(log(areaG),30);
n=N/sum(N);
% frequancy based on number of the grains
figure
bar(a,n);
title('Grain Area, N/sum(N)')
% frequency based on the area of the grains
figure
na=N.*exp(a)/sum(N.*exp(a));
bar(a,na)
title('Grain Area, n*exp(a)')

% image of the grain structure
figure
imshow(phi)
title(comstumstring)
figure
imshow(L)



% ************** Time Histrory Grains Statistics **************

% Read all the saved structures from the savedir folder at specified tn and
% then extract Area, Diameter and Perimeter of grains and put them in a
% cell array of A, D and P respectively. It takes long time depending on
% number of the saved frames
clear
savedir='set7';
load(strcat(pwd,'/',savedir,'/','setings','.mat'))
% figure
ni=0;
i=1;
for tn=[50:100 105:5:1000 1020:20:4000 4050:50:10000]
    filename=strcat(pwd,'/',savedir,'/',num2str(tn),'.mat');
    load(filename)
    phi(logical(ppf))=0;
    [areaG,diamG,perimG]=grainstat(phi);
    A{i}=areaG;
    D{i}=diamG;
    P{i}=perimG;
    i=i+1
end
% save all statistics
save(strcat(pwd,'/',savedir,'/','Grains','.mat'))

% ****** Plotings ******
% read statistics and find the average of the values, one can calculates
% other statistical properties inside the loop
clear
savedir='set8';
load(strcat(pwd,'/',savedir,'/','Grains','.mat'))
i=1;
tn=[50:100 105:5:1000 1020:20:4000 4050:50:10000];
tnV=tn(1:392);
for tn=tnV
    An(i)=mean(A{i});
    Dn(i)=sqrt(mean((D{i}).^2));
%   Dn(i)=mean(D{i});
    Pn(i)=mean(P{i});
    GBArea(i)=sum(P{i})/2;
    i=i+1;
end
tn=[50:100 105:5:1000 1020:20:4000 4050:50:10000];
tn=tnV;
plot(tn(1:end),An(1:end))
plot(tn(1:end),Pn(1:end))
plot(tn(1:end),GBArea(1:end))


cutn=50;
tn=tn(cutn:end);Dn=Dn(cutn:end);Pn=Pn(cutn:end);An=An(cutn:end);GBArea=GBArea(cutn:end);

% plotting Dn and fit
figure
% tn=[50:100 105:5:1000 1020:20:4000 4050:50:10000];
plot(tn,(Dn*delx).^2)
plot(tn(1:end),Dn(1:end)*delx)
pp=polyfit(tn,(Dn*delx).^2,1);
Dfit=sqrt(polyval(pp,tn));
hold on
plot(tn,Dfit,'r')
text((tn(end)+tn(1))/2,Dfit(end),strcat('k= ', num2str(pp(1))));
xlabel('time (time unit)')
ylabel('Average Diameter of Grains')

% plotting Dn and fit [Logs]
figure
% tn=[50:100 105:5:1000 1020:20:4000 4050:50:10000];
plot(log(tn),log(Dn*delx))
pp=polyfit(log(tn),log(Dn*delx),1);
Dfitlog=polyval([pp],log(tn));
hold on
plot(log(tn),Dfitlog,'r')
text(log((tn(end)+tn(1))/2),Dfitlog(end),strcat('n= ', num2str(pp(1))));
xlabel('Log time (time unit)')
ylabel('Log Average Diameter of Grains')


% plotting An and fit
figure
% tn=[50:100 105:5:1000 1020:20:4000 4050:50:10000];
plot(tn(1:end),An(1:end)*delx^2)
ppA=polyfit(tn,(An*delx^2),1);
Afit=polyval(ppA,tn);
hold on
plot(tn,Afit,'r')
text((tn(end)+tn(1))/2,Afit(end),strcat('k=',num2str(ppA(1))));
xlabel('time (time unit)')
ylabel('Average Area of Grains')

% plotting Pn and fit
figure
tn=[50:100 105:5:1000 1020:20:4000 4050:50:10000];
plot(tn(1:end),Pn(1:end)*delx)
pp=polyfit(tn,(Pn*delx).^2,1);
Pfit=sqrt(polyval(pp,tn));
hold on
plot(tn,Pfit,'r')
text((tn(end)+tn(1))/2,Pfit(end),strcat('k= ', num2str(pp(1))));


% energy dissipation
figure
sigma=1;
% plot(tn(1:end),1./(GBArea(1:end)*delx*sigma).^1.5)
Egb=(GBArea(1:end)*delx*sigma)/nboxsize/mboxsize
plot(tn(1:end),1./Egb.^2)
ppE=polyfit(tn,1./Egb.^2,1);
Efit=polyval(ppE,tn);   
hold on
plot(tn,Efit,'r')
text((tn(end)+tn(1))/2,Efit(end),strcat('k= ', num2str(ppE(1))));
xlabel('time (time unit)')
ylabel('1/Energy of the System^2')

