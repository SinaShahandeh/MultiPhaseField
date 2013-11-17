function [cur1,cur2,detach,cent1,cent2]=analysecontourgrains(eta,xparticle,yparticle,tn,ppf)
% For making white regions brighter
global mboxsize nboxsize
global delx
p=2;
% phi=sum(eta(:,:,1:p).^2,3);
% phi(phi>1)=nan;
% graymap=255/(max(max(phi-min(min(phi)))))*(phi-min(min(phi)));
% mag=1;
% subplot(3,1,1)
% % imshow(imresize(uint8(graymap),mag));
% imshow(uint8(graymap));
% hold on
% plot(xparticle*mag,yparticle*mag,'r.')
% title(strcat('Time= ', num2str(tn)))
% 
% contour(eta(:,:,1),[0.5 0.5],'color','b');
% contour(eta(:,:,2),[0.5 0.5],'color','r');
% 

% curvatures

% finding points where value of eta==0.5
x=[0:nboxsize-1]*delx;
y=[0:mboxsize-1]*delx;
% [X,Y]=meshgrid(x,y);
C1=contourc(x,y,eta(:,:,1),[0.5 0.5]);
C2=contourc(x,y,eta(:,:,2),[0.5 0.5]);

ind1=find(C1(1,:)==0.5);
ind2=find(C2(1,:)==0.5);

detach=0;
if size(ind2,2)==2
    detach=1;
end
% ind holds positions of points where two curve separation happens
if size(ind1,2)==1
    ind1(2)=size(C1,2);
else
    ind1(2)=ind1(2)-1;
end
if size(ind2,2)==1
    ind2(2)=size(C2,2);
else
    ind2(2)=ind2(2)-1;
end
% position of contours
xmid=nboxsize*delx/2;

x1=C1(1,2:ind1(2));
y1=C1(2,2:ind1(2));
x2=C2(1,2:ind2(2));
y2=C2(2,2:ind2(2));

% x1=x1(x1>50 & x1<150);
% y1=y1(x1>50 & x1<150);
% x2=x2(x2>50 & x2<150);
% y2=y2(x2>50 & x2<150);

% Polynomial fitting

l=length(x1);
x1=x1(fix(l/3):fix(2/3*l));
y1=y1(fix(l/3):fix(2/3*l));
l=length(x2);
x2=x2(fix(l/3):fix(2/3*l));
y2=y2(fix(l/3):fix(2/3*l));
try
    sp1=polyfit(x1,y1,4);
    spd1=polyder(sp1);
    cur1=polyval(polyder(spd1),xmid)/abs((1+polyval(spd1,xmid)^2)^1.5);
catch
    cur1=nan;
end
try
    sp2=polyfit(x2,y2,4);
    spd2=polyder(sp2);
    cur2=polyval(polyder(spd2),xmid)/abs((1+polyval(spd2,xmid)^2)^1.5);
catch
    cur2=nan;
end

% Spline Fitting
%
% try
%     sp1=spline(x1,y1);
%     spd1=splineder(sp1);
%     cur1=ppval(splineder(spd1),xmid)/abs((1+ppval(spd1,xmid)^2)^1.5);
% catch
%     cur1=nan;
% end
% try
%     sp2=spline(x2,y2);
%     spd2=splineder(sp2);
%     cur2=ppval(splineder(spd2),xmid)/abs((1+ppval(spd2,xmid)^2)^1.5);
% catch
%     cur2=nan;
% end

% h=subplot(3,1,3);
% curves=get(h,'children');
% % for the first value of the plot:
%
% if tn==1
%     tni=tn;
%     curi1=cur1;
%     curi2=cur2;
% else
%     tni=get(curves(1),'xData');
%     curi1=get(curves(2),'yData');
%     tni=[tni tn];
%     curi1=[curi1 cur1];
%     curi2=get(curves(1),'yData');
%     curi2=[curi2 cur2];0
% end
% hold off
%
% plot(tni,curi1,'b')
% hold on
% plot(tni,curi2,'r')

% SPEED
try
    cent1=interp1(x1,y1,xmid,'cubic');
catch
    cent1=nan;
end
try
    cent2=interp1(x2,y2,xmid,'cubic');
catch
    cent2=nan;
end

% 
% 
% axis equal
% axis([0 nboxsize 0 mboxsize])
% hold off
%
% subplot(3,1,2)
% xpos=fix(nboxsize/2);
% eta1=eta(:,xpos,1);
% eta2=eta(:,xpos,2);
% x=[1:mboxsize]*delx;
% plot(x,eta1,'b')
% hold on
% plot(x,eta2,'r')
% plot(x,ppf(:,xpos),'g')
% hold off

