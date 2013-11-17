function [Calculation_Error,averageitter]=implicit3(delt)
% clear all
% phase field parameters
L=1;
m=2;
gamma=1.5;
kappa=4;

G=[0 0.05];
% geometry settings
p=2;

delx=1;      % length unit per pixel
nboxsize=100/delx; % x axis in pixels

% delt=2;
timesteps=100/delt;
x1=20*delx;
eta=zeros(2,nboxsize);
xi=[0:nboxsize-1]*delx;
eta(1,:)=0.5*(1+tanh(sqrt(m/2/kappa)*(x1-xi)));
eta(2,:)=0.5*(1-tanh(sqrt(m/2/kappa)*(x1-xi)));

%% implicit
eta2=eta;
ittr=0;

for tn=1:timesteps
    eta2=eta; % future time step is similar to now
    eta2past=eta; error=100;
    itr=0;
    while (error>1e-7)
        %        eta2(1,1) = eta2(1,nboxsize-1);
        %        eta2(1,nboxsize) = eta2(1,2);
        %        eta2(2,1) = eta2(2,nboxsize-1);
        %        eta2(2,nboxsize) = eta2(2,2);
        %
        %        eta(1,1) = eta(1,nboxsize-1);
        %        eta(1,nboxsize) = eta(1,2);
        %        eta(2,1) = eta(2,nboxsize-1);
        %        eta(2,nboxsize) = eta(2,2);
        %
        eta2past=eta2;
        for p=1:2
            for j=2:nboxsize-1
                sumterm=eta2(1,j)^2+eta2(2,j)^2;
                ceta=(eta2(p,j));
                del2=1/delx^2*(eta2(p,j+1)+eta2(p,j-1));
                sumtermp=(sumterm-ceta.^2);
                pre=1/(1+L*kappa*2*delt/delx^2-L*delt*m+L*delt*6*G(p)+L*delt*m*2*gamma*sumtermp);
                eta2(p,j)=pre*(eta(p,j)-delt*L.*(m*(ceta^3)-kappa.*del2-6*ceta^2*G(p)));
            end
        end
        error=sum(sum((eta2-eta2past).^2));
        ittr=ittr+1;
        itr=itr+1;

%         plot(xi,eta(1,:))
%         hold on
%         plot(xi,eta(2,:),'r')
%         title(num2str(tn))
%         hold off
%         pause(0.01)
    end
    eta=eta2;
end
% figure
% plot(xi,eta(1,:),'.')
% hold on
% plot(xi,eta(2,:),'r')
% title(num2str(tn))
% hold off
% pause(0.01)

pos=interp1(eta(1,:),xi,0.5);
vel=(pos-x1)/(fix(timesteps)*delt);
mobility=3/2*L*sqrt(2*kappa/m);
analytical_vel=mobility*G(2);
Calculation_Error=(vel-analytical_vel)/analytical_vel*100;
averageitter=ittr/timesteps;

