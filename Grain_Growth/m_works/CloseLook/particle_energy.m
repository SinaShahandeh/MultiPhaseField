% calculating energy of the system on a circle arround the particle with 
% radius equalto diameter

function Epart=particle_energy(E,xparticle,yparticle,diameter,mboxsize,nboxsize)

%making a big circle with 1 with the same center as particle
ppf=zeros(mboxsize,nboxsize);

if mod(diameter,2)==0
    y1=yparticle-fix(diameter/2);y2=yparticle+fix(diameter/2)-1;
    x1=xparticle-fix(diameter/2);x2=xparticle+fix(diameter/2)-1;
else
    y1=yparticle-fix(diameter/2);y2=yparticle+fix(diameter/2);
    x1=xparticle-fix(diameter/2);x2=xparticle+fix(diameter/2);
end
if y1<1
    y1=1;
    disp('particle out of box')
end
if y2>nboxsize
    y2=nboxsize;
    disp('particle out of box')
end
if x1<1
    x1=1;
    disp('particle out of box')
end
if x2>mboxsize
    x2=mboxsize;
    disp('particle out of box')
end

part=imcircle(diameter);
onepixel=find(part==1);
ppf(y1:y2,x1:x2)=part(1:length(y1:y2),1:length(x1:x2));
% multiply ppf by E to choose only E values inside the diameter.

Epart=sum(sum(ppf.*E))/length(onepixel);
