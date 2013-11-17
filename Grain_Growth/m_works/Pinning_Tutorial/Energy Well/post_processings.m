%% post processings

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



pf=sum(sum(ppf))/mboxsize/nboxsize;
t=timevec(2:end);
figure
plot(t,eqpos,'r')


pp=polyfit(t,eqpos,1);
intvel=pp(1)


% ----- N particle_interface --------

% for loading the result of simulationi has been done before
% load('/media/disk/sim_res/Nparticle_interface/4/4.mat')
% particles volume fraction
pf=sum(sum(ppf))/mboxsize/nboxsize;
MV=MetaVol+pf.*MetaVol;
t=timevec(2:end);
figure
plot(t,MV,'r')


pp=polyfit(t,MV,1);
intvel=pp(1)/(nboxsize*delx)

intvel0=2.7725/2;
Pf=DelG(2)*(1-intvel/intvel0)


% 