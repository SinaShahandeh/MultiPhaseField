function drawgrains(phi,xparticle,yparticle,tn)
% For making white regions brighter
phi(phi>1)=nan;

graymap=255/(max(max(phi-min(min(phi)))))*(phi-min(min(phi)));
%         subplot(2,1,1)

mag=1;
if mag~=1
    imshow(imresize(uint8(graymap),mag));
else
    imshow(uint8(graymap));
end
title(strcat('Time= ', num2str(tn)))
hold on
plot(xparticle*mag,yparticle*mag,'r.')
hold off

