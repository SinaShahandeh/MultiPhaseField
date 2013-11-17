%%
clear
syms x R L y
R=solve('R-sqrt(R^2-L^2/16)-x',R)
gama=atan(y/(R-x))
beta=acos(y/(R^2+y^2));
alpha=

c=(R-x)^2+R^2-2*R*(R-x)*cos(alpha)
c=simple(c);
xprime=simple(sqrt(c^2-y^2))
xi=1
L=1
xprimei=subs(xprime,'L',L);
xprimei=subs(xprimei,'x',xi)
ezplot(xprimei,[-L/2 L/2])

