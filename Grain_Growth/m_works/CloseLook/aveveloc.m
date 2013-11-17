% take average of velocity field in a given area comes from phi and also
% calculates diameter of circular grain

function [vel]=aveveloc(Mveloc,phi,ppf,Mdetadt,shape)

global delx nboxsize mboxsize scale

%removing ppf from phi
phi(ppf==1)=1;

% % choosing only top of the dome
% if strcmp(shape,'dome')
% phi(:,[1:fix(nboxsize/2-15*scale/5) fix(nboxsize/2+15*scale/5):end])=1;
% end

%% selection based on phi:

% vel=mean(Mveloc(find(phi<0.7)));

%% selection based on detadt:
% choosing places where detadt is high
Mdetadt(ppf==1)=0;
maxM=0.9*max(max(abs(Mdetadt)));
vel=mean(Mveloc(find(abs(Mdetadt)>maxM)));




%% Finging the diameter of a circular grain from the phi field
% im=imfill(phi<0.65,'holes');
% L=bwlabel(im);
% D=regionprops(L,'EquivDiameter');
% D=[D.EquivDiameter]*delx;
