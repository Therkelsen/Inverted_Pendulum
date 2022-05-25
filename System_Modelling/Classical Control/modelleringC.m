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

cart = (mc+mp)*ddx + mp*l*ddtheta*cos(theta)-mp*l*sin(theta)*dtheta^2== F - bc*dx
pend = mp*l*ddtheta+mp*l*ddx*cos(theta)-sin(theta)*mp*l*dx*dtheta-mp*l*sin(theta)*(g+dx*dtheta) == -bp*dtheta

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
xSize = 750; ySize = 650;
xLeft = 100; yTop = 0;
set(fig,'Position',[xLeft yTop xSize ySize])
legend('Pendulum','Cart')

fig = figure()
stepplot(GP, GC)
title('System Step Responses')
xSize = 750; ySize = 650;
xLeft = 100; yTop = 0;
set(fig,'Position',[xLeft yTop xSize ySize])
xlabel('$Time~$','interpreter','latex')
ylabel('$Amplitude~$(meters)','interpreter','latex')
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
% xSize = 750; ySize = 650;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])
% 
% fig = figure()
% rlocus(GC*KsC)
% title('Cart Root Locus')
% xSize = 750; ySize = 650;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])

% kpP = 1;
% kiP = 1;
% kdP = 1;
% 
% kpC = 1;
% kiC = 1;
% kdC = 1;
% 
% KsP = kpP +kiP*(1/s) + kdP*s
% KsC = kpC +kiC*(1/s) + kdC*s
% 
% fig = figure()
% rlocus(GP*KsP)
% title('Pendulum Root Locus')
% xSize = 750; ySize = 650;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])
% 
% fig = figure()
% rlocus(GC*KsC)
% title('Cart Root Locus')
% xSize = 750; ySize = 650;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])

% kpP = -1;
% kiP = -1;
% kdP = -1;

kpC = -11.28;
kiC = -0.69;
kdC = 23.96;

% KsP = kpP +kiP*(1/s) + kdP*s
KsC = kpC +kiC*(1/s) + kdC*s

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

t1 = (0:0.001:0.35)';

% fig = figure()
% stepplot(HP,t1)
% title('PID Pendulum Step Response')
% xSize = 750; ySize = 650;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])
% xlabel('$Time~$','interpreter','latex')
% ylabel('$Amplitude~$(meters)','interpreter','latex')
% legend('Pendulum');
% grid on

LC = KsC*GC
HC = (LC)/(1+LC)

t2 = (0:0.001:0.1)';

% fig = figure()
% impulseplot(HP,t1)
% title('PID Pendulum Impulse Response')
% xSize = 750; ySize = 650;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])
% xlabel('$Time~$','interpreter','latex')
% ylabel('$Amplitude~$(meters)','interpreter','latex')
% legend('Pendulum');
% grid on

fig = figure()
stepplot(HP,HC,t1)
title('Parallel PID System Step Responses')
xSize = 750; ySize = 650;
xLeft = 100; yTop = 0;
set(fig,'Position',[xLeft yTop xSize ySize])
xlabel('$Time~$','interpreter','latex')
ylabel('$Amplitude~$(meters)','interpreter','latex')
legend('Pendulum','Cart');
grid on

% fig = figure()
% impulseplot(HP,HC,t2)
% title('Parallel PID System Impulse Responses')
% xSize = 750; ySize = 650;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])
% xlabel('$Time~$','interpreter','latex')
% ylabel('$Amplitude~$(meters)','interpreter','latex')
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

%% Modern Control
close all

A
B
C
D

states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'u'};
outputs = {'x'; 'phi'};

% This assumes full state feedback, which will not be the case irl
Q = C'*C;
Q(1,1) = 5000;
Q(3,3) = 100
R = 1;
K = lqr(A,B,Q,R)

Ac = [(A-B*K)];
Bc = [B];
Cc = [C];
Dc = [D];

states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'r'};
outputs = {'x'; 'phi'};

Cn = [1 0 0 0];
sys_ss = ss(A,B,Cn,0);
Nbar = rscale(sys_ss,K)

sys_cl = ss(Ac,Bc*Nbar,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

init_cond = [1; 2; 3; 4];
t = 0:0.01:5;
cart_ref = 0;
r = cart_ref * ones(size(t));
figure()
[y,t,x] = lsim(sys_cl, r, t, init_cond);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','$Cart~Position~$(m)','interpreter','latex')
set(get(AX(2),'Ylabel'),'String','$Pendulum~Angle~$(rad)','interpreter','latex')
xlabel('$Time~$(seconds)','interpreter','latex')
title('Step Response with Precompensation and LQR Control')

% Observer-based control is implemented as a replacement

ob = obsv(sys_ss);
observability = rank(ob)

poles = eig(Ac)

P = [-40 -41 -42 -43];
L = place(A',C',P)'

Ace = [(A-B*K) (B*K);
       zeros(size(A)) (A-L*C)];
Bce = [B*Nbar;
       zeros(size(B))];
Cce = [Cc zeros(size(Cc))];
Dce = [0;0];

states = {'x' 'x_dot' 'phi' 'phi_dot' 'e1' 'e2' 'e3' 'e4'};
inputs = {'r'};
outputs = {'x'; 'phi'};

sys_est_cl = ss(Ace,Bce,Cce,Dce,'statename',states,'inputname',inputs,'outputname',outputs);

init_cond = [init_cond; zeros(4,1)]
t = 0:0.01:5;
r = cart_ref*ones(size(t));
figure()
[y,t,x]=lsim(sys_est_cl, r, t, init_cond);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','$Cart~Position~$(m)','interpreter','latex')
set(get(AX(2),'Ylabel'),'String','$Pendulum~Angle~$(rad)','interpreter','latex')
xlabel('$Time~$(seconds)','interpreter','latex')

title('Step Response with Observer-Based State-Feedback Control')

    
function[Nbar] = rscale(a,b,c,d,k) 
    
    % Given the single-input linear system: 
    % . 
    % x = Ax + Bu 
    % y = Cx + Du
    % and the feedback matrix K,
    %
    % the function rscale(sys,K) or rscale(A,B,C,D,K) 
    % finds the scale factor N which will 
    % eliminate the steady-state error to a step reference 
    % for a continuous-time, single-input system 
    % with full-state feedback using the schematic below: 
    %
    %                     /---------\ 
    % R         +       u | .       | 
    % ---> N --->() ----> | X=Ax+Bu |--> y = Cx ---> y 
    %           -|        \---------/ 
    %            |             | 
    %            |<---- K <----| 
    % 
    % 8/21/96 Yanjie Sun of the University of Michigan 
    % under the supervision of Prof. D. Tilbury 
    % 6/12/98 John Yook, Dawn Tilbury revised error(nargchk(2,5,nargin));
    
    % --- Determine which syntax is being used --- nargin1 = nargin;
    
    if (nargin==2), % System form 
        
        [A,B,C,D] = ssdata(a); K=b;
    
    elseif (nargin==5), % A,B,C,D matrices A=a;
        A = a;
        B = b;
        C = c;
        D = d;
        K = k;
    else
        error('Input must be of the form (sys,K) or (A,B,C,D,K)')
    end;
    
    % compute Nbar s = size(A,1);
    s = size(A,1);
    Z = [zeros([1,s]) 1];
    N = inv([A,B;C,D])*Z';
    Nx = N(1:s);
    Nu = N(1+s);
    Nbar = Nu + K*Nx; 
end