clc, clear, close all
sympref('AbbreviateOutput', false);
syms theta(t) x(t)  u m M L g bcart bpend t I

  %HER LAVES MODELERING VHA. LAGRANGE 

Ekin = 0.5*(M+m)*diff(x, t)^2 + m*L*diff(x, t) * diff(theta, t)*cos(theta)+0.5*(m*L^2 + I)*diff(theta, t)^2
Epot = -m*g*L*cos(theta)

Lagrange = Ekin - Epot


l1 = diff(diff(Lagrange, diff(x, t)), t) - diff(Lagrange, x) == - bcart*diff(x, t) + u
l2 = diff(diff(Lagrange, diff(theta, t)), t) - diff(Lagrange, theta) == - bpend*diff(theta, t)

ddx = isolate(l2, diff(x, t, t))
ddt = isolate(l2, diff(theta, t, t))

la = subs(l1, diff(x, t, t), rhs(ddx));
eq2 = rhs(isolate(la, diff(theta, t, t)))

lb = subs(l1, diff(theta, t, t), rhs(ddt));
eq4 = rhs(isolate(lb, diff(x, t, t)))


eq1 = diff(x, t)
%eq2 = (m*L*sin(theta)*diff(theta, t)^2 - m*g*cos(theta)*sin(theta) - bcart*diff(x, t))/(M+m*sin(theta)^2) + u/(M+m*sin(theta)^2)
eq3 = diff(theta, t)
%eq4 =  ((M+m)*g*sin(theta) + bcart*cos(theta)*diff(x, t) - m*L*sin(theta)*cos(theta)*diff(theta, t)^2)/(L*(M+m*sin(theta)^2))  - (cos(theta)*u)/(L*(M+m*sin(theta)^2))                             % -(L*m*cos(theta)*diff(x, t) + L*g*m*sin(theta))/(I + L^2*m*cos(theta)^2 + L^2*m*sin(theta)^2)  % == omega'


% HER OPRETTES STATE-SPACE-VEKTORERNE

A_J = [diff(eq1, x), diff(eq1, diff(x, t)), diff(eq1, theta), diff(eq1, diff(theta, t));
    diff(eq2, x), diff(eq2, diff(x, t)), diff(eq2, theta), diff(eq2, diff(theta, t));
    diff(eq3, x), diff(eq3, diff(x, t)), diff(eq3, theta), diff(eq3, diff(theta, t));
    diff(eq4, x), diff(eq4, diff(x, t)), diff(eq4, theta), diff(eq4, diff(theta, t))]

B_J = [diff(eq1, u); diff(eq2, u); diff(eq3, u); diff(eq4, u)]


%HER UDFØRES LINEARISERING OMKRING PUNKTET (theta, diff(theta, t), diff(x, t)) = (pi, 0, 0)

A = subs(A_J, {x, theta, diff(theta, t), diff(x, t), u}, {0, pi, 0, 0, 0}) % pi = linearisering omkring opretstående punkt
B = subs(B_J, {x, theta, diff(theta, t), diff(x, t), u}, {0, pi, 0, 0, 0})% 0 =linearisering omkring nedadstående punkt
C = [1 0 0 0; 0 0 1 0] % vi måler x og theta (output) 
D = [0; 0]



% HER ER NOGLE AF DE FYSISKE STØRRELSER FOR SYSTEMET

 mPendul = 0.084; % masse af pendul [kg]
 varM =   0.5; % masse af vogn [kg]
 mStang = 0.082;  % masse af stang [kg]
 varm =  mStang + mPendul;
 varL =  0.35; % total længde af stang [m]
 varg = 9.82; % tyngdeaccelerationen [m/s^2]
% lConveyor = 1.72; %længde af conveyorbælte [m]
% rollerRadius = 0.05; %radius af bælte-rullerne [m]
 varbcart = 5; %dæmpning af conveyorbælte [N/(m/s)]
 varbpend = 0.0012; %dæmpning for pendul [N/(m/s)]


 % IERTIMOMENTET FOR PENDULET

 varI = (1.0/3.0)*varm*varL^2;


%DER INDSÆTTES NUMERISKE VÆRDIER HER

A = double(subs(A, {m, M, L, g, bcart, bpend, I}, {varm, varM, varL, varg, varbcart, varbpend, varI}))
B = double(subs(B, {m, M, L, g, bcart, bpend, I}, {varm, varM, varL, varg, varbcart, varbpend, varI}))
C = double(subs(C, {m, M, L, g, bcart, bpend, I}, {varm, varM, varL, varg, varbcart, varbpend, varI}))
D = double(subs(D, {m, M, L, g, bcart, bpend, I}, {varm, varM, varL, varg, varbcart, varbpend, varI}))


s = tf('s');
tf = C*(inv(s*eye(4, 4)-A))*B+D


%% HER OPRETTES OVERFØRINGSFUNKTIONERNE

clear


tfC = tf([3.957 78.24 26.69 1538 -119.5 7551  0], [1 39.58 386.4  283.1  7698  -1036 37750  0  0])
tfP = tf([1.847 36.55 -4.924 358.8 0], [1 39.58 386.4  283.1  7698  -1036   37750])

t = 0:0.01:10;
%u = zeros(size(t));
%u(t == 1) = 1;

u = ones(size(t));

plot(u)

figure()
lsim(tfP, u, t)
title("Pendulum simulation")

figure()
lsim(tfC, u, t)
title("Cart simulation")

