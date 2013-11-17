%% time derivative and spline interpolation
i=12; % y axis
j=22; % x axis
t=timevec(2:end);
eta35=reshape(MMeta1(i,j,1:end),1,length(t));

fn=spaps(t,eta35,1e-11);
detadt=fnder(fn);
% time -time derivative plot
figure
fnplt(detadt)
hold on
plot(t,reshape(MMdetadt(i,j,1:end),1,length(t)),'r--')
title('Time derivative of the eta function at i and j node point')
% time function plot
figure
fnplt(fn)
hold on
plot(t,reshape(MMeta1(i,j,1:end),1,length(t)),'r--')
title('The eta fucntion and its spline fit')





x=1:44;
y=1:40;
Fn=csapi({x,y},eta1');   % cubic
% Fn=spaps({x,y},eta1',1e-9);  % smoothing


%% position derivative
[nablaetax,nablaetay]=gradient(MMeta1(:,:,100),delx,delx);

x=1:44;
y=1:40;
Fn=csapi({x,y},eta1');   % cubic
Fn=spaps({x,y},eta1',1e-1);  % smoothing

Fnd=fndir(Fn,eye(2));

xi=linspace(1,44,200);
yi=linspace(1,40,200);
[X,Y]=meshgrid(xi,yi);

Fneta = reshape(fnval(Fn,[X(:) Y(:)].'), length(yi),length(xi));
Fngrads = reshape(fnval(Fnd,[X(:) Y(:)].'),2,length(yi),length(xi));
 
% mesh(X,Y,squeeze(Fngrads(1,:,:)))
% mesh(X,Y,squeeze(Fngrads(2,:,:)))

% plot of the function
figure
eta=MMeta1(:,:,100);
plot(eta(:,j),'linewidth',2)
hold on
plot(yi,Fneta(:,end/2),'r-.')

% plot of the derivative
figure
plot(nablaetay(:,j))
hold on
plot(yi,squeeze(Fngrads(2,:,end/2)),'m')
title(' Y derivative')

figure
plot(nablaetax(:,j))
hold on
plot(yi,squeeze(Fngrads(1,:,end/2)),'m')
title(' X derivative')

%% Velocity Field





%% Velocity from iso-surface
t=timevec(2:end);
centpp=spaps(t,Mcent1,1e-3);
Velpp=fnder(centpp);

figure
fnplt(centpp)
hold on
plot(t,Mcent1,'r-.')
title('Position of Interface from iso-\eta')

figure
fnplt(Velpp)
hold on
plot(t(2:end),diff(Mcent1)/delt,'r-.')
title('Velocity from iso-\eta')

%% Volume Change
Mvol=MetaVol(2,:)/(mboxsize*nboxsize*delx^2);
t=timevec(2:end);
Volpp=spaps(t,Mvol,1e-7);
Velpp=fnder(Volpp);

figure
fnplt(Volpp)
hold on
plot(t,Mvol,'r-.')
title('Volume of the phase field')

figure
fnplt(Velpp)
hold on
plot(t(2:end),diff(Mvol)/delt,'r-.')
title('dV/dt')

