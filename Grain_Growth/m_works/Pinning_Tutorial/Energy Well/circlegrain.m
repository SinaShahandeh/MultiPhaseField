function [ppf] =circlegrain...
    (mboxsize,nboxsize,xp,yp,diameter,shape)
% xp and yp is position of the circle_centre
ppf=zeros(mboxsize,nboxsize);

% y1,x1 and y2,x2 are bounding box of the circle
if mod(diameter,2)==0
    y1=yp-fix(diameter/2);y2=yp+fix(diameter/2)-1;
    x1=xp-fix(diameter/2);x2=xp+fix(diameter/2)-1;
else
    y1=yp-fix(diameter/2);y2=yp+fix(diameter/2);
    x1=xp-fix(diameter/2);x2=xp+fix(diameter/2);
end
x1p=x1; % previous x1
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
    x2=nboxsize;
end

xcut=x1-x1p;
part=imcircle(diameter);
% part(part==0)=-1;
ppf(y1:y2,x1:x2)=part(1:length(y1:y2),xcut+1:xcut+length(x1:x2));

% to make in a dome shape
if strcmp(shape,'dome')
    ppf(fix((y1+y2)/2):end,x1:x2)=1;
end
