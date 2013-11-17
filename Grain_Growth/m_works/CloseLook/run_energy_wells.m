
%% ---- ENERGY WELL ---- 
%% LOOP for energy of the system with one particle and interface
% % 2. A grain boundary with particle in the middle
% clear all
% global nboxsize mboxsize
% global delx
% scale=4;
% epsilon=linspace(1,10,5)
% for ieps=1:20
%     x=linspace(0.2,0.8,20);
%     for ni=1:length(x)
%         [eta,ppf,ME,E]=stationarymin_p2(x(ni),epsilon(ieps));
%         MEM(:,:,ni)=ME;
%         MinterfaceE(ni)=E/(nboxsize*delx);
%     end
%     WellProfile(ieps,:)=MinterfaceE;
%     % measuring well depth
%     % interpolation of well
%     Wellx=linspace(0,mboxsize*delx,1000);
%     Welly=interp1(fix(mboxsize*x)*delx,WellProfile,Wellx,'cubic');
%     WellDepth(ieps)=Minterface(1)-min(Welly);
%     
% end
% 
% figure
% plot(fix(mboxsize*x)*delx,WellProfile)

%% looping for the 2 independent variables of epsilon and particle size
clear all
global nboxsize mboxsize
global delx
scale=3;
epsilon=linspace(1e-18,10e-18,5)
Pr=linspace(5,1,5);
for ipr=1:length(Pr)
    for ieps=1:length(epsilon)
        x=[0.2 linspace(0.3,0.7,15) 0.8];
        for ni=1:length(x)
            [eta,ppf,ME,E]=stationarymin_p2(x(ni),epsilon(ieps),Pr(ipr));
            MEM(:,:,ni)=ME;
            MinterfaceE(ni)=E;
        end
        WellProfile(ieps,ipr,:)=MinterfaceE;
        % measuring well depth
        % interpolation of well
        Wellx=linspace(0,mboxsize*delx,1000);
        Welly=interp1(fix(mboxsize*x)*delx,MinterfaceE,Wellx,'cubic');   
        WellDepth(ieps,ipr)=MinterfaceE(1)-min(Welly);
        lowerindex=find(Welly<(MinterfaceE(1)-0.05*WellDepth(ieps,ipr)));
        WellWidth(ieps,ipr)=Wellx(lowerindex(end))-Wellx(lowerindex(1));
        ipr
        ieps
    end
    save(['Potential_Well_results/ipr' num2str(ipr) '.mat'])
end
beep



