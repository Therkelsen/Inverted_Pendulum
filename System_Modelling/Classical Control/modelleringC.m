%% System modelling
clc, clear, close all
sympref('AbbreviateOutput', false);
syms mc mp l theta dtheta ddtheta x dx ddx F bc bp Jp g;

g_ = 9.82;    %Earth gravitational force

mc_ = 0.5;    %Mass of cart
bc_ = 5;      %Friction coeficient of cart
%bc_ = 0;

mp_ = 0.084;  %Mass of pendulum
bp_ = 0.0012; %Friction coeficcient of pendulum
%bp_ = 0;
l_ = 0.35;    %length of pendulum

%Jp_ = (1/3)*mp_*l_^2
Jp_ = 0;

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
sys_tf = tf(sys_ss);

% Individual transfer functions
[num, denum] = tfdata(sys_tf);

s = tf('s');

GC = tf(num(1),denum(1))
GP = tf(num(2),denum(2));

fig = figure()
pzmap(GP, GC)
title('System PZ Map')
xSize = 900; ySize = 800;
xLeft = 100; yTop = 0;
set(fig,'Position',[xLeft yTop xSize ySize])
legend('Pendulum','Cart')

fig = figure()
stepplot(GP, GC)
title('System Step Responses')
xSize = 900; ySize = 800;
xLeft = 100; yTop = 0;
set(fig,'Position',[xLeft yTop xSize ySize])
xlabel('$Time~[s]$','interpreter','latex')
ylabel('$Amplitude~[m]$','interpreter','latex')
legend('Pendulum', 'Cart');
grid on

%% Manual Parallel PID tuning using root locus
s = tf('s');

% kpP = 1;
% kiP = 0;
% kdP = 0;
% 
% kpC = 1;
% kiC = 0;
% kdC = 0;
% 
% KsP = kpP +kiP*(1/s) + kdP*s
% KsC = kpC +kiC*(1/s) + kdC*s
% 
% fig = figure()
% rlocus(GP*KsP)
% title('Pendulum Root Locus')
% xSize = 900; ySize = 800;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])
% 
% fig = figure()
% rlocus(GC*KsC)
% title('Cart Root Locus')
% xSize = 900; ySize = 800;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])

kpP = 1;
kiP = 1;
kdP = 1;

kpC = 1;
kiC = 1;
kdC = 1;

KsP = kpP +kiP*(1/s) + kdP*s
KsC = kpC +kiC*(1/s) + kdC*s

fig = figure()
rlocus(GP*KsP)
title('Pendulum Root Locus')
xSize = 900; ySize = 800;
xLeft = 100; yTop = 0;
set(fig,'Position',[xLeft yTop xSize ySize])

fig = figure()
rlocus(GC*KsC)
title('Cart Root Locus')
xSize = 900; ySize = 800;
xLeft = 100; yTop = 0;
set(fig,'Position',[xLeft yTop xSize ySize])

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

%% Auto tuning of Parallel PID controllers
%kp = proportional gain, ki = integral grain
% kd = derivative gain, Tf = filter order
[KsP,infoP] = pidtune(GP, 'PIDF')
[kpP,kiP,kdP,TfP] = piddata(KsP);

[KsC,infoC] = pidtune(GC, 'PIDF')
[kpC,kiC,kdC,TfC] = piddata(KsC);

%% Control structure analysis
close all
LP = KsP*GP
HP = (LP)/(1+LP)

t1 = (0:0.001:2)';

fig = figure()
stepplot(HP,t1)
title('PID Pendulum Step Response')
xSize = 900; ySize = 800;
xLeft = 100; yTop = 0;
set(fig,'Position',[xLeft yTop xSize ySize])
xlabel('$Time~[s]$','interpreter','latex')
ylabel('$Amplitude~[m]$','interpreter','latex')
legend('Pendulum');
grid on

LC = KsC*GC
HC = (LC)/(1+LC)

t2 = (0:0.001:0.1)';

% fig = figure()
% impulseplot(HP,t1)
% title('PID Pendulum Impulse Response')
% xSize = 900; ySize = 800;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])
% xlabel('$Time~[s]$','interpreter','latex')
% ylabel('$Amplitude~[m]$','interpreter','latex')
% legend('Pendulum');
% grid on

% fig = figure()
% stepplot(HP,HC,t1)
% title('Parallel PID System Step Responses')
% xSize = 900; ySize = 800;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])
% xlabel('$Time~[s]$','interpreter','latex')
% ylabel('$Amplitude~[m]$','interpreter','latex')
% legend('Pendulum','Cart');
% grid on

% fig = figure()
% impulseplot(HP,HC,t2)
% title('Parallel PID System Impulse Responses')
% xSize = 900; ySize = 800;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])
% xlabel('$Time~[s]$','interpreter','latex')
% ylabel('$Amplitude~[m]$','interpreter','latex')
% legend('Pendulum','Cart');
% grid on

%% Digitization of controllers
Ts = 0.001;

z = tf('z');

sReplace = (2/Ts)*((z-1)/(z+1));

KzP = kpP + kiP*(1/sReplace) + kdP*sReplace

KzC = kpC + kiC*(1/sReplace) + kdC*sReplace

%% Cascade controllers

kpP = 695;
kiP = 3.8e+03;
kdP = 31.8;

KsP = kpP +kiP*(1/s) + kdP*s

kpC = 0.00468134655401907;
kiC = 0.000103385823399229;
kdC = 0.0505001089730261;
TfC = 1.41006414988528;

KsC = kpC + kiC*(1/s) + kdC*(TfC/(1+TfC*(1/s)))

%%
close all

A
B
C
D


Obs = [C; C*A; C*A*A; C*A*A*A]

Con = [B  A*B   A*A*B   A*A*A*B]

nul = [0, 0; 0, 0; 0, 0; 0, 0]
nul1 = [0, 0; 0, 0]
nul2 = [0; 0]

Ae = [A, nul;C, nul1]
Be = [B; nul2]
Ce = [C, nul1]

rank(Obs)

rank(Con)

p = [-1.1 -1.2 -100 -110 -120 -130];
K = place(Ae,-Be,p)

Fi = K(:, [5, 6])
F = K(:, 1:4)

nul3 = [0, 0, 0]

p1 = [-50, -40, -45, -47]
L = place(A', -C', p1)
L = L'


