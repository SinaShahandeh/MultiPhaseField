% distribution could be 'normal' or 'uniform'
function [ppf,xparticle,yparticle] =particledistro...
    (nboxsize,mboxsize,particles_number,diameter,distribution)

ppf=zeros(nboxsize,mboxsize);
if particles_number==0
    xparticle=nan;
    yparticle=nan;
end

% first determine particles position based on spacial distribution


%% for chunk of normal distributed particles in non-uniform distribution
if strcmp(distribution,'normal')
chunknum=60;
xparticle=[];yparticle=[];
for chunki=1:chunknum
    xchunk=mboxsize*rand(1,1);
    ychunk=nboxsize*rand(1,1);
    x=fix(30*randn(1,fix(particles_number/chunknum))+xchunk)+1;
    y=fix(30*randn(1,fix(particles_number/chunknum))+ychunk)+1;
    xparticle=[xparticle x];
    yparticle=[yparticle y];
end
% imposing periodic boundary condition on particles outside the domain
xparticle(xparticle>mboxsize)=mboxsize-xparticle(xparticle>mboxsize);
xparticle(xparticle<1)=mboxsize+xparticle(xparticle<1);
yparticle(yparticle>nboxsize)=nboxsize-yparticle(yparticle>nboxsize);
yparticle(yparticle<1)=nboxsize+yparticle(yparticle<1);
end

%% For uniform distribution of particles
if strcmp(distribution,'uniform')
yparticle=fix(nboxsize*rand(1,particles_number))+1;
xparticle=fix(mboxsize*rand(1,particles_number))+1;
end 

% now put particles into ppf
for np=1:particles_number
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
end


