% save desired variables to disk with matlab workspace format

function  savegrains(phi,xparticle,yparticle,tn,costumstring)

filename=strcat(pwd,'/',costumstring,'/',num2str(tn),'.mat')
save(filename)

