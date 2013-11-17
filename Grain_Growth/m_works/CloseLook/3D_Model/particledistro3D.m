function [ppf,xp,yp,zp] =particledistro3D...
    (nboxsize,mboxsize,lboxsize,particles_number,diameter)

ppf=zeros(mboxsize,nboxsize,lboxsize);
if particles_number==0
    xp=nan;
    yp=nan;
    zp=nan;
end

for np=1:particles_number
    xp=fix(nboxsize/2);
    yp=fix(mboxsize/2*9/10);
    zp=fix(lboxsize/2);
if mod(diameter,2)==0
    y1=yp-fix(diameter/2);y2=yp+fix(diameter/2)-1;
    x1=xp-fix(diameter/2);x2=xp+fix(diameter/2)-1;
    z1=zp-fix(diameter/2);z2=zp+fix(diameter/2)-1;
else
    y1=yp-fix(diameter/2);y2=yp+fix(diameter/2);
    x1=xp-fix(diameter/2);x2=xp+fix(diameter/2);
    z1=zp-fix(diameter/2);z2=zp+fix(diameter/2);
end
if y1<1
    y1=1;
end
if y2>mboxsize
    y2=mboxsize;
end
if x1<1
    x1=1;
end
if x2>nboxsize
    x2=naboxsize;
end
if z1<1
    z1=1;
end
if z2>lboxsize
    z2=laboxsize;
end

    part=imsphere(diameter/2);
    ppf(y1:y2,x1:x2,z1:z2)=part(1:length(y1:y2),1:length(x1:x2),1:length(z1:z2));
end








% for np=1:particles_number
%     yparticle(np)=fix(nboxsize*rand(1,1))+1;xparticle(np)=fix(mboxsize*rand(1,1))+1;
%     circle=imcircle(diameter);
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

