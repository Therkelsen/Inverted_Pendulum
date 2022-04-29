clc, clear, close all
sympref('AbbreviateOutput', false);
syms theta(t) x(t)  u m M l g bcart bpend t I


  %HER lAVES MODElERING VHA. lAGRANGE 

Ekin = 0.5*(M+m)*diff(x, t)^2 + m*l*diff(x, t) * diff(theta, t)*cos(theta)+0.5*(m*l^2 + I)*diff(theta, t)^2
Epot = -m*g*l*cos(theta)

lagrange = Ekin - Epot

l1 = (M+m)*diff(x, t, t) + m*l*diff(theta, t, t) == u - bcart*diff(x, t)
l2 = (m*l^2 + I)*diff(theta, t, t) + m*l*diff(x, t, t) + m*g*l*theta == 0


clc, clear, close all

 mPendul = 0.084; % masse af pendul [kg]
 M =   0.5; % masse af vogn [kg]
 mStang = 0.082;  % masse af stang [kg]
 m =  mStang + mPendul;
 l =  0.35; % total længde af stang [m]
 g = -9.82; % tyngdeaccelerationen [m/s^2]
 b = 5; %dæmpning af conveyorbælte [N/(m/s)]

  I = (1.0/3.0)*m*l^2;

s =tf('s')

q = (M+m)*(M*l + I) - (m*l)^2

tfP = (-m*l*s/q)/(s^3 + ((b*(m*l^2 + I))/q)*s^2 + ((M+m)*(g*m*l)/q)*s + (b*m*g*l)/q)
tfC = ((m*l^2 + I)*s^2 + g*m*l/q)/(s^4 + (b*(m*l + I)/q)*s^3 + ((M+m)*g*m*l/q)*s^2 + (b*m*g*l/q)*s)

[nump, denp] = tfdata(tfP)
[numc, denc] = tfdata(tfC)

b1 = -(m*l^2 + I)*b/q
b2 = m^2*g*l^2/q
b3 = m*l*b/q
b4 = -(M+m)*m*g*l/q
w1 = (I+m*l^2)/q
w2 = -m*l/q

A = [0, 1, 0, 0; 0, b1, b2, 0; 0, 0, 0, 1; 0, b3, b4, 0]
B = [0; w1; 0; w2]
C = [1, 0, 0, 0; 0, 0, 1, 0]
D = [0; 0]

sys = ss(A, B, C, D)

t = 0:0.1:10
u = ones(size(t))
x0 = [0; 0; 0; 0]
% 
% figure()
% lsim(sys, u, t, x0)
% 
% 
% figure()
% impulse(tfP)
% 
% figure()
% impulse(tfC)

%[pidP info] = pidtune(tfP, 'pidf')

%pidP = -1.37*1000 -4.33*1000/s -93.2*(s/0.00198*s + 1);

%pidTuner(tfP, 'PIDF')

Ti = 0.1;
Td = 0.1;

pidP = -(1 + (1/Ti)*(1/s) + Td*s)

figure(1)
rlocus(tfP*pidP)

Ti = 2;
Td = 0.1;


K = 572*pidP

H = K*tfP/(1 + K*tfP)

figure(2)
lsim(H, u, t)

figure(3)
margin(K*tfP)






