%% explicit method convergence

Mdelt=logspace(-3,-2.5,10);
for n=1:length(Mdelt)
    MError(n)=explicit(Mdelt(n));
end
figure
semilogx(Mdelt,MError)
    


%% implicit method convergence

Mdelt=logspace(log10(1e-1),log10(10),100);
for n=1:length(Mdelt)
    [MError(n), MItr(n)]=implicit1(Mdelt(n));
end
figure
semilogx(Mdelt,MError)
ylabel('Error (%)')

% figure
% semilogx(Mdelt,MItr)
% ylabel('Number of Itterations')

figure
gain=Mdelt/0.05./MItr;
semilogx(Mdelt,gain)
ylabel('Speed up Gain')

%% implicit method convergence

Mdelt=logspace(log10(1e-1),log10(10),100);
for n=1:length(Mdelt)
    [MError(n), MItr(n)]=implicit2(Mdelt(n));
end
figure
semilogx(Mdelt,MError)
ylabel('Error (%)')

% figure
% semilogx(Mdelt,MItr)
% ylabel('Number of Itterations')

figure
gain=Mdelt/0.05./MItr;
semilogx(Mdelt,gain)
ylabel('Speed up Gain')

%% implicit method convergence

Mdelt=logspace(log10(1e-1),log10(10),100);
for n=1:length(Mdelt)
    [MError(n), MItr(n)]=implicit3(Mdelt(n));
end
figure
semilogx(Mdelt,MError)
ylabel('Error (%)')

% figure
% semilogx(Mdelt,MItr)
% ylabel('Number of Itterations')

figure
gain=Mdelt/0.05./MItr;
semilogx(Mdelt,gain)
ylabel('Speed up Gain')

%% implicit method convergence

Mdelt=logspace(log10(1e-1),log10(10),100);
for n=1:length(Mdelt)
    [MError(n), MItr(n)]=implicit4(Mdelt(n));
end
figure
semilogx(Mdelt,MError)
ylabel('Error (%)')

% figure
% semilogx(Mdelt,MItr)
% ylabel('Number of Itterations')

figure
gain=Mdelt/0.05./MItr;
semilogx(Mdelt,gain)
ylabel('Speed up Gain')

%% implicit method convergence

Mdelt=logspace(log10(1e-1),log10(10),100);
for n=1:length(Mdelt)
    [MError(n), MItr(n)]=implicit5(Mdelt(n));
end
figure
semilogx(Mdelt,MError)
ylabel('Error (%)')

% figure
% semilogx(Mdelt,MItr)
% ylabel('Number of Itterations')

figure
gain=Mdelt/0.05./MItr;
semilogx(Mdelt,gain)
ylabel('Speed up Gain')

%% implicit method convergence

Mdelt=logspace(log10(1e-1),log10(10),100);
for n=1:length(Mdelt)
    [MError(n), MItr(n)]=implicit6(Mdelt(n));
end
figure
semilogx(Mdelt,MError)
ylabel('Error (%)')

% figure
% semilogx(Mdelt,MItr)
% ylabel('Number of Itterations')

figure
gain=Mdelt/0.05./MItr;
semilogx(Mdelt,gain)
ylabel('Speed up Gain')

