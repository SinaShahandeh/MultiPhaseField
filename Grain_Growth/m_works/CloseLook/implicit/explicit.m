function [out]=explicit(delti)
% phase field parameters
L=1;
m=2;
gamma=1.5;
kappa=4;
epsilon=5;

G=[0 0.05];
% setings structure
settings.L=L(1);
settings.alpha=m;
settings.beta=m;
settings.gamma=gamma(1);
settings.kappa=kappa(1);
settings.epsilon=epsilon(1);
settings.DelG=G;
settings.accuracy='low';
% geometry settings
p=2;
delx=1;      % length unit per pixel
nboxsize=100/delx; % x axis in pixels


 delt=0.1;
% delt=delti
timesteps=100/delt;
x1=20*delx;
eta=zeros(2,nboxsize);
xi=[0:nboxsize-1]*delx;
eta(1,:)=0.5*(1+tanh(sqrt(m/2/kappa)*(x1-xi)));
eta(2,:)=0.5*(1-tanh(sqrt(m/2/kappa)*(x1-xi)));
%% explicit forward euler

eta2=eta;
for tn=1:timesteps
    for j=2:nboxsize-1
        sumterm=eta(1,j)^2+eta(2,j)^2;
        for p=1:2

            del2=1/delx^2*(eta(p,j+1)+eta(p,j-1)-2*eta(p,j));
            sumtermp=eta(p,j)*sumterm-eta(p,j).^3;
            detadt=-L.*(m*(-eta(p,j)+eta(p,j)^3+2*gamma*sumtermp)-kappa.*del2+6*eta(p,j)*(1-eta(p,j))*G(p));
            eta2(p,j)=eta(p,j)+delt*detadt;
        end
    end
eta=eta2;

plot(eta(1,:))
hold on
plot(eta(2,:),'r')
title(num2str(tn))
hold off
pause(0.05)
end
figure
plot(xi,eta(1,:))
hold on
plot(xi,eta(2,:),'r')
title(num2str(tn))
hold off
pause(0.01)

pos=interp1(eta(1,:),xi,0.5);
vel=(pos-x1)/(timesteps*delt);
mobility=3/2*L*sqrt(2*kappa/m);
analytical_vel=mobility*G(2);
Calculation_Error=(vel-analytical_vel)/analytical_vel*100;

out=Calculation_Error;
%% implicit backward euler
% figure
% eta=zeros(2,nboxsize);
% eta2=eta;
% eta(1,1:x)=1;
% eta(2,x+1:end)=1;
% ittr=0;
% for tn=1:timestepss
%     eta2=eta; % future time step is similar to now
%     eta2past=eta; error=100;
%     while error >0.000001
%         eta2past=eta2;
%         for j=2:nboxsize-1
%             sumterm=eta(1,j)^2+eta(2,j)^2;
%             for p=1:2
%                 ceta=eta2(p,j);
%                 del2=1/delx^2*(eta2(p,j+1)+eta2(p,j-1)-2*ceta);
%                 sumtermp=ceta*sumterm-ceta.^3;
%                 detadt=-L.*(m*(-ceta+ceta^3+2*gamma*sumtermp)-kappa.*del2+6*ceta*(1-ceta)*G(p));
%                 eta2(p,j)=eta(p,j)+delt*detadt;
%             end
%         end
%         error=sum(sum((eta2-eta2past).^2));
%         ittr=ittr+1;
%     end
%     eta=eta2;
% end
% toc
% 
% averageitter=ittr/timesteps
% plot(eta(1,:))
% hold on
% plot(eta(2,:),'r')
% title(num2str(tn))
% hold off
% pause(0.01)
