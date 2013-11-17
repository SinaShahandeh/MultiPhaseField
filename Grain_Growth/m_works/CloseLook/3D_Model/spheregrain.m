function [ggf] =spheregrain...
    (mboxsize,nboxsize,lboxsize,xp,yp,zp,diameter)
% xp and yp is position of the circle_centre
ggf=zeros(mboxsize,nboxsize,lboxsize)-1;


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

part=imsphere(fix(diameter/2));
part(part==0)=-1;
ggf(y1:y2,x1:x2,z1:z2)=part(1:length(y1:y2),1:length(x1:x2),1:length(z1:z2));
% to make in a dome shape
for i=1:lboxsize
    [yi,xi]=find(ggf(:,:,i)==1);
    xi1=min(xi);xi2=max(xi);
    ggf(fix(yp):end,xi1:xi2,i)=1;
end
