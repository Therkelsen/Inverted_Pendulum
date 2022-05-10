%%  HER lAVES MODElERING VHA. lAGRANGE 

clc, clear, close all
sympref('AbbreviateOutput', false);
syms theta(t) x(t)  u m M l g bcart bpend t I


Ekin = 0.5*(M+m)*diff(x, t)^2 + m*l*diff(x, t) * diff(theta, t)*cos(theta)+0.5*(m*l^2 + I)*diff(theta, t)^2
Epot = -m*g*l*cos(theta)

lagrange = Ekin - Epot

l1 = diff(diff(lagrange, diff(x, t)), t) - diff(lagrange, x) == - bcart*diff(x, t) + u
l2 = diff(diff(lagrange, diff(theta, t)), t) - diff(lagrange, theta) == 0


l1 = (M+m)*diff(x, t, t) + m*l*diff(theta, t, t) == u - bcart*diff(x, t)
l2 = (m*l^2 + I)*diff(theta, t, t) + m*l*diff(x, t, t) + m*g*l*theta == 0

%% state feedback

clc, clear, close all

 mPendul = 0.084; % masse af pendul [kg]
 M =   0.5; % masse af vogn [kg]
 mStang = 0.082;  % masse af stang [kg]
 m =  mStang + mPendul;
 l =  0.35; % to tal længde af stang [m]
 g = 9.82; % tyngdeaccelerationen [m/s^2]
 b = 5; %dæmpning af conveyorbælte [N/(m/s)]

  I = (1.0/3.0)*m*l^2;

s =tf('s')

q = (M+m)*(M*l + I) - (m*l)^2

tfP = tf([1, 0, 0], [m*l-((m+M)/(m*l))*(m*l^2+ I), 0, -g*(m+M) - b*(m*l^2 + I)/m*l, -b*g])
tfC = tf([-(m*l^2+I)/(m*l), 0, -g], [m*l-(m+M)*(m*l^2+I)/(m*l), -((m*l^2+I)*b)/(m*l), g*(m+M), -g*b, 0])

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

t = 0:0.1:20
u = ones(size(t))

% x0 = [0; 0; 0; 0]
% 
% figure()
% lsim(sys, u, t, x0)

kp = -1;
ki = -1;
kd = -1;

K = kp +ki*(1/s) + kd*s

figure()
rlocus(tfP*K)
title('Pendulum Root Locus')

[K,info] = pidtune(tfP, 'PIDF')

[kp,ki,kd,Tf] = piddata(K)

Tp = (K*tfP)/(1+K*tfP);
Tc = (K*tfC)/(1+K*tfC);

impulseplot(Tp, Tc)
legend('Pendulum response','Cart response')

% kp = 100;
% ki = 1000;
% kd = 50;
% 
% K = kp + ki*(1/s) + kd*s

%Tp = feedback(K*tfP,1)

% figure()
% impulse(Tp)
% title('Pendulum Impulse Response')
% figure()
% bode(Tp)
% title('Pendulum Bode Plot')

% Tc = feedback(K*tfC,1)
% figure()
% impulse(Tc)
% title('Cart Impulse Response')
% figure()
% bode(Tc)
% title('Cart Bode Plot')

% pole(tfC)
% pole(tfP)
% 
% figure()
% rlocus(tfP*K)
% title('Pendulum Root Locus')
% 
% figure()
% rlocus(tfC*K)
% title('Cart Root Locus')

LC = tfC*K

% K = 120*K;

Ts = 0.001;

z = tf('z');

sReplace = (2/Ts)*((z-1)/(z+1));

discreteTek = kp + ki*(1/sReplace) + kd*sReplace

%Hp = tfP/(1+tfP*K)

%Hc = tfC/(1+tfP*K)

% figure()
% impulse(Hp, t)
% title('Pendulum Impulse Response')
% 
% figure()
% impulse(Hc, t)
% title('Cart Impulse Response')
% 
% figure()
% margin(K*tfP)
% title('Pendulum Stability Margins')
% 
% figure()
% margin(K*tfC)
% title('Cart Stability Margins')

