%grain growth model: Fan and Chen, acta mater 1997, Vol 45, p 611
%Fourier spectral method: Trefethen, L.N., Spectral Methods in MATLAB,
%SIAM, Society for Industrial and Applied Mathematics, 2000, Philadelphia,
%PA
%Fourier spectral method: Chen, L.Q. and Shen, J., Comp. Phys. Comm., 1998,
%108, 147


%solution PDE: fft in space, semi-implicit Euler method in time
%N = number of gridpoints in each D
%k1 = L*deltat (=0.25 in fan chen)
%k2 = kappa/deltax*deltax (= 0.5 in fan chen)
%order_par = number of order parameters (36 in fan chen)
%n = number of timesteps
%pl = 'plot' or 'save'
%pl_times = vector with times to save
%n_old: simulation starts at time step n_old + 1; if n_old = -1  ([sumetasqu,eta] = grgr3d(N,k1,k2,m,order_par,pl,pl_times,n,-1);) => simulation starts from random structure
N=[50 50 10];
k1=0.5;
k2=0.4;
m=2;
order_par=36;
n=1000;
pl='plot'
pl_times=[1:1000];
n_old=-1;

% function [sumetasqu,eta] = grgr3d(N,k1,k2,m,order_par,pl,pl_times,n,n_old,eta);

j=1; %used in loop for plotting, saving
%construction vectors in fourier space
g = (2*pi/N(1)) * [0:1:N(1)/2 -N(1)/2+1:1:-1];
gg = (2*pi/N(2)) * [0:1:N(2)/2 -N(2)/2+1:1:-1];
ggg = (2*pi/N(3)) * [0:1:N(3)/2 -N(3)/2+1:1:-1];
[g1,g2,g3] = meshgrid(gg,g,ggg); %N*N*N-matrix
clear g gg ggg

%distributing nuclei in the matrix
nuclei_fraction=0.05;
if n_old == -1
    for i = 1 : order_par
        eta{i} =fix(rand(N(1),N(2),N(3))+nuclei_fraction);
    end 
    n_old = 0;
end

ETA = cell(1,order_par); sumetasqu = zeros(N(1),N(2),N(3));
for i = 1 : order_par
    sumetasqu = sumetasqu+eta{i}.^2; %sum over all squared phase fields
    ETA{i} = fftn(eta{i}); %fft of orderparameters
end

XB =  k1*k2*(g1.^2+g2.^2+g3.^2) + 1; 
bulk = complex(zeros(N(1),N(2),N(3))); %%fft of bulk driving force 

for k = n_old + 1:n
    tic
    if mod(k,50) == 0
    fprintf('iterationstep %g \n',k)
    end
    for i = 1 : order_par
        a = eta{i};
        b = a.^2;
        bulk = fftn(m*a.*( - b + 2*sumetasqu - 1)); %fft of bulk driving force (eta^3-eta+2*eta*sum_{j\neq i}(eta) is rearranged a little to reduce number of computations)
        %semi-implicit euler for time discretisation 
        ETA{i}=(ETA{i}-k1*bulk)./XB; %fft of orderparameters at time n+1
        eta{i} = real(ifftn(ETA{i})); % orderparameter at time n+1
        
        eta{i}(eta{i}>1)=1;
        eta{i}(eta{i}<0)=0;
    end
    
    sumetasqu = eta{1}.^2; 
    for i = 2 : order_par
        sumetasqu = sumetasqu + eta{i}.^2;
    end

    if pl=='plot'
        if j <= length(pl_times) & k == pl_times(j)
            j = j + 1;
%             figure(1)
%             clf
            pcolor(sumetasqu(:,:,1));
            caxis([0 1]);
            colormap gray
            shading('interp')
            axis off
            box off 
            set(gca,'DataAspectratio',[1 1 1])
            pause(0.01);
        end
    end

    if pl=='save'
        if j <=length(pl_times) & k==pl_times(j)
           j=j+1;
             save(strcat(invoer, 'r', int2str(k)),'sumetasqu');
        end
    end
    
    
    toc
end


if pl =='save'
    save(strcat('grain_growth','einde',int2str(n)),'eta','sumetasqu')
end


