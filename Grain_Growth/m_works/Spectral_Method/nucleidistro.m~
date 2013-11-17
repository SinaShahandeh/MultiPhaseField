
function [etapf] =nucleidistro3d(nboxsize,mboxsize,kboxsize,particles_number,diameter,distribution)

etapf=zeros(nboxsize,mboxsize,kboxsize);


% first determine particles position based on spacial distribution

%% For uniform distribution of particles
if strcmp(distribution,'uniform')
yparticle=fix(nboxsize*rand(1,particles_number))+1;
xparticle=fix(mboxsize*rand(1,particles_number))+1;
zparticle=fix(kboxsize*rand(1,particles_number))+1;

end 
% cubic nuclei
etapf(xparticle,yparticle,zparticle)=1;
% 
% % now put particles into ppf
% for np=1:length(yparticle)
%     % checking if the particle is truncated with the border and then truncate
%     % it. It actually should move the truncated part to the other side of the
%     % priodic boundary condition but it does not. Most of the time particles
%     % are just one grid point ant this ifs are unusable.
%     if mod(diameter,2)==0
%         y1=yparticle(np)-fix(diameter/2);y2=yparticle(np)+fix(diameter/2)-1;
%         x1=xparticle(np)-fix(diameter/2);x2=xparticle(np)+fix(diameter/2)-1;
%     else
%         y1=yparticle(np)-fix(diameter/2);y2=yparticle(np)+fix(diameter/2);
%         x1=xparticle(np)-fix(diameter/2);x2=xparticle(np)+fix(diameter/2);
%     end
%     if y1<1
%         y1=1;
%     end
%     if y2>nboxsize
%         y2=nboxsize;
%     end
%     if x1<1
%         x1=1;
%     end
%     if x2>mboxsize
%         x2=mboxsize;
%     end
%     part=imcircle(diameter);
%     etapf(y1:y2,x1:x2)=part(1:length(y1:y2),1:length(x1:x2));
% end
% 
