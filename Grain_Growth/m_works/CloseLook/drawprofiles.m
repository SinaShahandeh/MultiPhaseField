function drawprofiles(eta,ppf)
% For making white regions brighter
global mboxsize nboxsize
global delx

xpos=fix(nboxsize/2);
eta1=eta(:,xpos,1);
eta2=eta(:,xpos,2);
x=[1:mboxsize]*delx;
plot(x,eta1,'b')
hold on
plot(x,eta2,'r')
plot(x,ppf(:,xpos),'g')
hold off