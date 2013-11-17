%This function generates a phere with radious r in the containig box
function [sph]=imsphere(r)

sph=zeros(2*r,2*r,2*r);
[x,y,z]=sphere(r*10);
x=fix(x*2*r/2)+r+1;
y=fix(y*2*r/2)+r+1;
z=fix(z*2*r/2)+r+1;
for i=1:size(x,1)
    for j=1:size(x,2);
    sph(x(j,i),y(j,i),z(j,i))=1;
    end
end

sph=imfill(sph);


