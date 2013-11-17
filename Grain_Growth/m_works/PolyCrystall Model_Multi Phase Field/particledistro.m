function [ppf,xparticle,yparticle] =particledistro...
    (nboxsize,mboxsize,particles_number,diameter)

ppf=zeros(nboxsize,mboxsize);
if particles_number==0
    xparticle=nan;
    yparticle=nan;
end

%% Random Distribution of Particles
% for np=1:particles_number
%     yparticle(np)=fix(nboxsize*rand(1,1))+1;
%     xparticle(np)=fix(mboxsize*rand(1,1))+1;
%     circle=imcircle(diameter);
%  % checking if the particle is truncated with the border and then truncate
%  % it. It actually should move the truncated part to the other side of the
%  % priodic boundary condition but it does not. Most of the time particles
%  % are just one grid point ant this ifs are unusable.
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
%     ppf(y1:y2,x1:x2)=part(1:length(y1:y2),1:length(x1:x2));
% end

%% Arranged quare distribution of particles
np=1;
rownum=fix(sqrt(particles_number));
columnnum=fix(sqrt(particles_number));
 for npi=0:rownum
    for npj= 0: columnnum
    yparticle(np)=fix(npi/rownum*nboxsize)+1;
    xparticle(np)=fix(npj/columnnum*mboxsize)+1;
    circle=imcircle(diameter);
    
 % checking if the particle is truncated with the border and then truncate
 % it. It actually should move the truncated part to the other side of the
 % priodic boundary condition but it does not. Most of the time particles
 % are just one grid point ant this ifs are unusable.
    if mod(diameter,2)==0
        y1=yparticle(np)-fix(diameter/2);y2=yparticle(np)+fix(diameter/2)-1;
        x1=xparticle(np)-fix(diameter/2);x2=xparticle(np)+fix(diameter/2)-1;
    else
        y1=yparticle(np)-fix(diameter/2);y2=yparticle(np)+fix(diameter/2);
        x1=xparticle(np)-fix(diameter/2);x2=xparticle(np)+fix(diameter/2);
    end
    if y1<1
        y1=1;
    end
    if y2>nboxsize
        y2=nboxsize;
    end
    if x1<1
        x1=1;
    end
    if x2>mboxsize
        x2=mboxsize;
    end
    part=imcircle(diameter);
    ppf(y1:y2,x1:x2)=part(1:length(y1:y2),1:length(x1:x2));
    np=np+1;
    end
end

