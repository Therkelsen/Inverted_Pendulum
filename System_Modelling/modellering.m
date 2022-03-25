clc, clear, close all
sympref('AbbreviateOutput', false);
syms theta(t) x(t)  u m M L g bcart bpend t I


% HER LAVES MODELERING VHA. LAGRANGE D'ALEMBERT

Ekin = 0.5*(M*diff(x, t)^2+m*(x+L*diff(theta, t)*cos(theta))^2+m*(L*diff(theta, t)*sin(theta))^2+I*diff(theta, t)^2)
Epot = -m*g*L*cos(theta)

Lagrange = Ekin - Epot

% BEVÆGELSESLIGNINGERNE

eq1 = diff(x, t)
eq2 = diff(diff(Lagrange, diff(x, t)), t)-diff(Lagrange, x) + bcart*diff(x, t) - u
eq3 = diff(theta, t)
eq4 = diff(diff(Lagrange, diff(theta, t)), t)-diff(Lagrange, theta) + bpend*diff(theta, t)

% HER OPRETTES STATE-SPACE-VEKTORERNE

A = [diff(eq1, x), diff(eq1, diff(x, t)), diff(eq1, theta), diff(eq1, diff(theta, t));
    diff(eq2, x), diff(eq2, diff(x, t)), diff(eq2, theta), diff(eq2, diff(theta, t));
    diff(eq3, x), diff(eq3, diff(x, t)), diff(eq3, theta), diff(eq3, diff(theta, t));
    diff(eq4, x), diff(eq4, diff(x, t)), diff(eq4, theta), diff(eq4, diff(theta, t))]

B = [diff(eq1, u); diff(eq2, u); diff(eq3, u); diff(eq4, u)]


%HER UDFØRES LINEARISERING OMKRING PUNKTET (theta, diff(theta, t), diff(x, t)) = (pi, 0, 0)

A = subs(A, {x, theta, diff(theta, t), diff(x, t), u}, {0, pi, 0, 0, 0}) % pi = linearisering omkring opretstående punkt
B = subs(B, {x, theta, diff(theta, t), diff(x, t), u}, {0, pi, 0, 0, 0}) % 0 =linearisering omkring nedadstående punkt
C = [1 0 0 0; 0 0 1 0] % vi måler x og theta (output) 
D = [0; 0]



% HER ER NOGLE AF DE FYSISKE STØRRELSER FOR SYSTEMET

 mPendul = 0.084; % masse af pendul [kg]
 varM = 0.5; % masse af vogn [kg]
 mStang = 0.082;  % masse af stang [kg]
 varm = mStang + mPendul;
 varL = 0.35; % total længde af stang [m]
 varg = -9.82; % tyngdeaccelerationen [m/s^2]
% lConveyor = 1.72; %længde af conveyorbælte [m]
% rollerRadius = 0.05; %radius af bælte-rullerne [m]
 varbcart = 5; %dæmpning af conveyorbælte [N/(m/s)]
 varbpend = 0.0012; %dæmpning for pendul [N/(m/s)]


 % IERTIMOMENTET FOR PENDULET

 varI = (1.0/3.0)*m*L^2;


%DER INDSÆTTES NUMERISKE VÆRDIER HER

A = double(subs(A, {m, M, L, g, bpend, bcart, I}, {varm, varM, varL, varg, varbpend, varbcart, varI}))
B = double(subs(B, {m, M, L, g, bpend, bcart, I}, {varm, varM, varL, varg, varbpend, varbcart, varI}))
C = double(subs(C, {m, M, L, g, bpend, bcart, I}, {varm, varM, varL, varg, varbpend, varbcart, varI}))
D = double(subs(D, {m, M, L, g, bpend, bcart, I}, {varm, varM, varL, varg, varbpend, varbcart, varI}))


%HER OPRETTES OVERFØRINGSFUNKTIONERNE
s = tf('s');
tf = C*(inv(s*eye(4, 4)-A))*B+D

