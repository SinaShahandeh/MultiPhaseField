function [finaleta]=boundarycond(eta,ppf,value)

%fixed BC on ppf points
finaleta=eta;
finaleta(find(ppf==1))=value;