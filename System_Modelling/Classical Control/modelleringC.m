%% System modelling
clc, clear, close all
sympref('AbbreviateOutput', false);
syms mc mp l theta dtheta ddtheta x dx ddx F bc bp Jp g;

g_ = 9.82;    %Earth gravitational force

mc_ = 0.5;    %Mass of cart
bc_ = 5;      %Friction coeficient of cart

mp_ = 0.084;  %Mass of pendulum
bp_ = 0.0012; %Friction coeficcient of pendulum
%bp_ = 0;
l_ = 0.35;    %length of pendulum
Jp_ = (1/3)*mp_*l_^2
%Jp_ = 0;

%Ours (probably wrong)
%cart = (mc+mp)*ddx + mp*l*cos(theta)*ddtheta - mp*l*sin(theta)*dtheta^2 == F - bc*dx
%pend = (mp*l + Jp)*ddtheta + mp*l*ddx*cos(theta)+mp*g*l*sin(theta) == 0

cart = (mc+mp)*ddx-mp*l*cos(theta)*ddtheta+mp*l*sin(theta)*dtheta^2 == F - bc*dx
pend = (mp*l + Jp)*ddtheta-mp*l*ddx*cos(theta)-mp*g*l*sin(theta) == -bp*dtheta

sol = solve([cart,pend],ddx,ddtheta);

x1 = dx
x2 = sol.ddx
x3 = dtheta
x4 = sol.ddtheta

%Differentiate into jacobian:
An = [diff(x1,x) diff(x1,dx) diff(x1,theta) diff(x1,dtheta);
      diff(x2,x) diff(x2,dx) diff(x2,theta) diff(x2,dtheta);
      diff(x3,x) diff(x3,dx) diff(x3,theta) diff(x3,dtheta);
      diff(x4,x) diff(x4,dx) diff(x4,theta) diff(x4,dtheta)];

Bn = [diff(x1,F);
      diff(x2,F);
      diff(x3,F);
      diff(x4,F)];

%Linearize about 0 degrees, assuming zero speed. 

Al = subs(An, {theta dtheta x dx F},{0 0 0 0 0})
Bl = subs(Bn, {theta dtheta x dx F},{0 0 0 0 0})

%Insert Values
A = double(subs(Al, {Jp mc mp l bc bp g}, {Jp_ mc_ mp_ l_ bc_ bp_ g_}))
B = double(subs(Bl, {Jp mc mp l bc bp g}, {Jp_ mc_ mp_ l_ bc_ bp_ g_}))

vpa(A);
vpa(B);

%Values for matrix [cart pos, cart velocity, pend pos, pend velocity]
%velocity
C = double([1 0 0 0; 0 0 1 0]) %Cart speed
D = double(0);


% Statespace object
sys_ss = ss(double(A),double(B),double(C),D);

% Transferfunction from statespace
sys_tf = tf(sys_ss)

% Individual transfer functions
[num, denum] = tfdata(sys_tf);

GC = tf(num(1),denum(1));
GP = tf(num(2),denum(2));

%pole(GP)
%title('Pendulum PZ Map')
%pole(GC)
%title('Cart PZ Map')

%% Manual PID tuning using root locus
% kpP = -1;
% kiP = -1;
% kdP = -1;
%
% kpC = -1;
% kiC = -1;
% kdC = -1;
% 
% KsP = kpP +kiP*(1/s) + kdP*s
% KsC = kpC +kiC*(1/s) + kdC*s
%
% figure()
% rlocus(GP*KsP)
% title('Pendulum Root Locus')
% 
% figure()
% rlocus(GC*KsC)
% title('Cart Root Locus')
%
% kpP = -1;
% kiP = -1;
% kdP = -1;
%
% kpC = -1;
% kiC = -1;
% kdC = -1;
% 
% KsP = kpP +kiP*(1/s) + kdP*s
% KsC = kpC +kiC*(1/s) + kdC*s

%% Auto tuning of PID controllers
%kp = proportional gain, ki = integral grain
% kd = derivative gain, Tf = filter order
[KsP,infoP] = pidtune(GP, 'PID')
[kpP,kiP,kdP,TfP] = piddata(KsP);

[KsC,infoC] = pidtune(GC, 'PID')
[kpC,kiC,kdC,TfC] = piddata(KsC);

%% Control structure analysis
LP = KsP*GP;
LC = KsC*GC;

HP = (LP)/(1+LP);
HC = (LC)/(1+LC);

t1 = (0:0.001:1)';
t2 = (0:0.001:0.05)';

figure()
stepplot(HP, 'b',t1)
title('Pendulum Step Response')

figure()
stepplot(HC,'r',t2)
title('Cart Step Response')

% figure()
% impulseplot(HP, 'b')
% title('Pendulum Impulse Response')
% 
% figure()
% impulseplot(HC,'r')
% title('Cart Impulse Response')

% figure()
% bode(Tp)
% title('Pendulum Bode Plot')
% 
% figure()
% bode(Tc)
% title('Cart Bode Plot')

% figure()
% margin(K*GP)
% title('Pendulum Stability Margins')
% 
% figure()
% margin(K*GC)
% title('Cart Stability Margins')

%% Digitization of controllers
Ts = 0.001;

z = tf('z');

sReplace = (2/Ts)*((z-1)/(z+1));

KzP = kpP + kiP*(1/sReplace) + kdP*sReplace

KzC = kpC + kiC*(1/sReplace) + kdC*sReplace
