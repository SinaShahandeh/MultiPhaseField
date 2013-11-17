% Find fraction of triple juinction that coinsice with particle position
function [tripleind,particleind]=TJonP(triples,particles)
% Mdistance is a matrix that calculates distance of each particle from
% tripple junctions
Mdistance=zeros(length(triples),length(particles));
for i=1:length(triples)
    for j=1:length(particles)
        dx=triples(1,i)-particles(1,j);
        dy=triples(2,i)-particles(2,j);
        Mdistance(i,j)=sqrt(dx^2+dy^2);
    end
end

coincidedist=2*sqrt(2);
% coincidedist=2*sqrt(3);
[tripleind,particleind]=find(Mdistance<coincidedist);
