function  savegrains(eta,xparticle,yparticle,zparticle,tn,savedir)

filename=strcat(pwd,'/',savedir,'/',num2str(tn),'.mat')
save(filename)
