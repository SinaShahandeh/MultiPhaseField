function [ppf,xparticle,yparticle] =particledistroN...
    (nboxsize,mboxsize,particles_fraction,diameter,ratio)

Sr=1/4;
ppf=zeros(nboxsize,mboxsize);
%generating particles number for different region
Hparticles_number=fix(particles_fraction/(pi*(diameter/2)^2)*nboxsize*Sr*mboxsize*ratio);
Lparticles_number=fix(particles_fraction/(pi*(diameter/2)^2)*nboxsize*(1-Sr)*mboxsize);
if particles_fraction==0
    xparticle=nan;
    yparticle=nan;
end

for np=1:Hparticles_number
    yparticle(np)=fix(nboxsize*Sr*rand(1,1)+((1-Sr)/2)*nboxsize)+1;
    xparticle(np)=fix(mboxsize*rand(1,1))+1;
    %putting particle in the middle of the xparticle and yparticle
    if mod(diameter,2)==0
        y1=yparticle(np)-fix(diameter/2);y2=yparticle(np)+fix(diameter/2)-1;
        x1=xparticle(np)-fix(diameter/2);x2=xparticle(np)+fix(diameter/2)-1;
    else
        y1=yparticle(np)-fix(diameter/2);y2=yparticle(np)+fix(diameter/2);
        x1=xparticle(np)-fix(diameter/2);x2=xparticle(np)+fix(diameter/2);
    end
    % check if particle is not out of the simulation box
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
    % making circular particle and putting it on the ppf
    part=imcircle(diameter);
    ppf(y1:y2,x1:x2)=part(1:length(y1:y2),1:length(x1:x2));
end

np=np+1;

while np<Lparticles_number+Hparticles_number+1
    yparticle(np)=fix(nboxsize*rand(1,1))+1;
    xparticle(np)=fix(mboxsize*rand(1,1))+1;
    if yparticle(np)<((1-Sr)/2)*nboxsize | yparticle(np)>((1-Sr)/2+Sr)*nboxsize
        %putting particle in middle of the xparticle and yparticle
        if mod(diameter,2)==0
            y1=yparticle(np)-fix(diameter/2);y2=yparticle(np)+fix(diameter/2)-1;
            x1=xparticle(np)-fix(diameter/2);x2=xparticle(np)+fix(diameter/2)-1;
        else
            y1=yparticle(np)-fix(diameter/2);y2=yparticle(np)+fix(diameter/2);
            x1=xparticle(np)-fix(diameter/2);x2=xparticle(np)+fix(diameter/2);
        end
        % check if particle is not out of the simulation box
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
        % making circular particle and putting it on the ppf
        part=imcircle(diameter);
        ppf(y1:y2,x1:x2)=part(1:length(y1:y2),1:length(x1:x2));
        np=np+1;

    end
end