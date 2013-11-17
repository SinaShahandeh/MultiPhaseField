function Khist=HistCurvature(savedir,m,kappa,Pz,tn,maxK,maxY,PlotType)

load ([savedir 'Khist_' num2str(tn) '.mat']);
Khist(abs(Khist)>0.05)=[];
Khist=abs(Khist);
Mtn_meanK=mean(abs(Khist));
intenergy=1/3*sqrt(2*m*kappa);



if strcmp(PlotType,'bar')==true;
    x=linspace(0,maxK+0.005,50); % bins
[N,X]=histc(Khist,x); %
    bar(x,N/length(Khist),'group');
    title([ 'Time step = ' num2str(tn) ], 'FontSize',24);
    h=gca;
    set(h,'LineWidth',2);
    set(h,'FontSize',20)
    xlabel('K','FontSize',24);
    ylabel('Length Fraction','FontSize',24);% pause(1)
    hold on
    
plot([Mtn_meanK Mtn_meanK], [0 maxY],'g-' ,'LineWidth',3)
plot(([Pz/intenergy Pz/intenergy]), [0 maxY],'r-' ,'LineWidth',2)
hold off

    set(gca,'Xlim',[0 maxK]);
set(gca,'Ylim',[0 maxY]);
end

if strcmp(PlotType,'plot')==true;
    x=linspace(0,maxK/Mtn_meanK,50); % bins
    [N,X]=histc(Khist/Mtn_meanK,x); %
    plot(x,N/length(Khist),'LineWidth',2);
    title([ 'Time step = ' num2str(tn) ], 'FontSize',16);
    h=gca;
    set(h,'LineWidth',1.5);
    set(h,'FontSize',14)
    xlabel('K/<K>','FontSize',16);
    ylabel('Length Fration','FontSize',16);% pause(1)
    hold on
    axis([min(x) maxK/Mtn_meanK 0 maxY])
end



