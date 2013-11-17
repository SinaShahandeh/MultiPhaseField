%% grain boundary gradient

L=1;
m=2;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);

syms x
eta=0.5*(1+tanh(sqrt(m/2/kappa)*x))
subplot(2,1,1)
ezplot(eta,[-5 5])
grid on

grad=diff(eta)
subplot(2,1,2)
ezplot(grad,[-5 5])
grid on

%% solute friction
figure
L=1;
m=2;
kappa=2;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
a=0.2;
b=4;
syms v;
Pf=a*v/(a+b*v^2);

% ezplot(Pf,[0 5])
vi=linspace(0,5,100);
Pfi=subs(Pf,vi);

plot(vi/mobility/intenergy,Pfi/intenergy)
ylabel('P_f/ \sigma_{gb}');
xlabel('V / (M \sigma_{gb})');
hold on

%% solving analytical curve for cahn drag
figure
clear
L=1;
m=2;
kappa=4;
mobility=3/2*L*sqrt(2*kappa/m);
intenergy=1/3*sqrt(2*m*kappa);
a=0.2;
b=4;
syms v;
Pf=a*v/(a+b*v^2);

delG=linspace(0,1,100);
vel_pure=mobility*delG;
plot(delG/intenergy,vel_pure/mobility/intenergy,'r')
hold on
n=0
for dG=delG
    n=n+1;
    e1=mobility*(dG-Pf)-v;
    e1sol=solve(e1);
    vel_sol(n)=eval(e1sol(1));
end
hold on
plot(delG/intenergy,vel_sol/mobility/intenergy,'b')
xlabel('\Delta G/ \sigma_{gb}');
ylabel('V / (M \sigma_{gb})');


