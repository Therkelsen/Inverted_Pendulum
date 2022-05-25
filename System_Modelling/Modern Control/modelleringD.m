%% System modelling
clc, clear, close all
sympref('AbbreviateOutput', false);
syms mc mp l theta dtheta ddtheta x dx ddx F bc bp Jp g t X(t) Theta(t)

g_ = 9.82    %Earth gravitational force

mc_ = 0.5    %Mass of cart
bc_ = 5      %Friction coeficient of cart

mp_ = 0.084  %Mass of pendulum
bp_ = 0.0012 %Friction coeficcient of pendulum
l_ = 0.35    %length of pendulum    (179 mm correct)

cart = (mc+mp)*ddx + mp*l*ddtheta*cos(theta)-mp*l*cos(theta)*dtheta^2== F - bc*dx
pend = mp*l^2*ddtheta + mp*l*ddx*cos(theta)-sin(theta)*mp*l*dx*dtheta - mp*l*sin(theta)*(g-dx*dtheta) == -bp*dtheta 

% lagrange = 0.5*(mc + mp)*diff(X, t)^2 + mp*diff(X, t)*l*diff(Theta, t)*cos(Theta) + 0.5*mp*l^2*diff(Theta, t)^2 + mp*g*l*cos(Theta)
% Wai = diff(diff(lagrange, diff(Theta, t)), t) - diff(lagrange, Theta)
% Wai1 = diff(diff(lagrange, diff(X, t)), t) - diff(lagrange, X)


sol = solve([cart,pend],ddx,ddtheta);

x1 = dx;
x2 = sol.ddx
x3 = dtheta;
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


%Linearize about 0 degrees, asuming zero speed. 

Al = subs(An, {theta dtheta x dx F},{pi 0 0 0 0})
Bl = subs(Bn, {theta dtheta x dx F},{pi 0 0 0 0})

%Insert Values
A = subs(Al, {mc mp l bc bp g}, {mc_ mp_ l_ bc_ bp_ g_});
B = subs(Bl, {mc mp l bc bp g}, {mc_ mp_ l_ bc_ bp_ g_});

vpa(A);
vpa(B);
A = double(A);
B = double(B);

%Values for matrix [cart pos, cart velocity, pend pos, pend velocity]
%velocity
C = [1 0 0 0; 0 0 1 0]; %Cart speed
D = 0;

% Statespace object
sys_ss = ss(double(A),double(B),double(C),D)

% Transferfunction from statespace
sys_tf = tf(sys_ss)

% Individual transfer functions
[num, denum] = tfdata(sys_tf);

GC = tf(num(1),denum(1))
GP = tf(num(2),denum(2))

% PLC sample time, used for digitization
Ts = 0.001;
% Saturation is decided by the motor max force
sat = 1168;

%Here I remove the numerical error on sys_tf_pend
%sys_tf_pend = tf([1.803, 0], [1, 9.058, 10.71, 88.53])

% fig = figure()
% pzmap(GP, GC)
% title('System PZ Map')
% xSize = 750; ySize = 650;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])
% legend('Pendulum','Cart')
% 
% fig = figure()
% stepplot(GP, GC)
% title('System Step Responses')
% xSize = 750; ySize = 650;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])
% xlabel('$Time~$','interpreter','latex')
% ylabel('$Amplitude~$(meters)','interpreter','latex')
% legend('Pendulum', 'Cart');
% grid on

%% Manual Parallel PID tuning using root locus
close all
s = tf('s');

kpP = 1;
kiP = -3460;
kdP = -9.07;

kpC = 1;
kiC = 1.9;
kdC = 0.12;

KsP = kpP +kiP*(1/s) + kdP*s
KsC = kpC +kiC*(1/s) + kdC*s

fig = figure()
rlocus(GP*KsP)
title('Pendulum Root Locus')
xSize = 750; ySize = 650;
xLeft = 100; yTop = 0;
set(fig,'Position',[xLeft yTop xSize ySize])

fig = figure()
rlocus(GC*KsC)
title('Cart Root Locus')
xSize = 750; ySize = 650;
xLeft = 100; yTop = 0;
set(fig,'Position',[xLeft yTop xSize ySize])

kpP = -354;
kiP = -3460;
kdP = -9.07;

kpC = 1;
kiC = 1.9;
kdC = 0.12;

KsP = kpP+kiP*(1/s)+kdP*s;
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

t1 = (0:0.001:1)';
t2 = (0:0.001:10)';

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

% fig = figure()
% stepplot(HC,t1)
% title('PID Cart Step Response')
% xSize = 750; ySize = 650;
% xLeft = 100; yTop = 0;
% set(fig,'Position',[xLeft yTop xSize ySize])
% xlabel('$Time~$','interpreter','latex')
% ylabel('$Amplitude~$(meters)','interpreter','latex')
% legend('Cart');
% grid on

fig = figure()
stepplot(HP,HC,t2)
title('Parallel PID System Step Responses')
xSize = 750; ySize = 650;
xLeft = 100; yTop = 0;
set(fig,'Position',[xLeft yTop xSize ySize])
xlabel('$Time~$','interpreter','latex')
ylabel('$Amplitude~$(rad \& meters)','interpreter','latex')
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

%% Modern Control
close all
Anew = [0,  1.0000, 0,  0; 
        -10.3406,   -0.0430,    0,  -9.0155;
        0, 0, 0, 1;
        -0.5206,    -0.0022,    0,  -9.0155]

A = Anew

% This assumes full state feedback, which might not be the case irl
% Without weights
Q = C'*C;
R = 1;
K = lqr(A,B,Q,R)

Ac = [(A-B*K)];
Bc = [B];
Cc = [C];
Dc = [D];

states = {'theta' 'theta_dot' 'x' 'x_dot'};
inputs = {'r'};
outputs = {'theta';'x'};

init_cond = [pi/8; 0; 1; 0];
t = 0:0.001:15;
cart_ref = 0;

sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

r = cart_ref * ones(size(t));
figure()
[y,t,x] = lsim(sys_cl, r, t, init_cond);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','$Cart~Position~$(meters)','interpreter','latex')
set(get(AX(2),'Ylabel'),'String','$Pendulum~Angle~$(radians)','interpreter','latex')
xlabel('$Time~$(seconds)','interpreter','latex')
title('Step Response with LQR Control')

% With precompensation and added weights in the Q matrix
cart_weight = 1500;
pend_weight = 5000;

Q = C'*C;
Q(1,1) = cart_weight
Q(3,3) = pend_weight
R = 1;
K = lqr(A,B,Q,R)

Ac = [(A-B*K)];
Bc = [B];
Cc = [C];
Dc = [D];


Cn = [0 0 1 0];
sys_ss = ss(A,B,Cn,0);
Nbar = rscale(sys_ss,K)

sys_cl = ss(Ac,Bc*Nbar,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

figure()
[y,t,x] = lsim(sys_cl, r, t, init_cond);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','$Cart~Position~$(meters)','interpreter','latex')
set(get(AX(2),'Ylabel'),'String','$Pendulum~Angle~$(radians)','interpreter','latex')
xlabel('$Time~$(seconds)','interpreter','latex')
title('Step Response with weighed LQR Control')

% Observer-based control is implemented as a replacement
% 
% ob = obsv(sys_ss);
% observability = rank(ob)
% 
% poles = eig(Ac)
% 
% P = [-40 -41 -42 -43];
% L = place(A',C',P)'
% 
% Ace = [(A-B*K) (B*K);
%        zeros(size(A)) (A-L*C)];
% Bce = [B*Nbar;
%        zeros(size(B))];
% Cce = [Cc zeros(size(Cc))];
% Dce = [0;0];
% 
% states = {'theta' 'theta_dot' 'x' 'x_dot' 'e1' 'e2' 'e3' 'e4'};
% inputs = {'r'};
% outputs = {'theta'; 'theta_dot'};
% 
% sys_est_cl = ss(Ace,Bce,Cce,Dce,'statename',states,'inputname',inputs,'outputname',outputs);
% 
% init_cond = [init_cond; zeros(4,1)]
% figure()
% [y,t,x]=lsim(sys_est_cl, r, t, init_cond);
% [AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
% set(get(AX(1),'Ylabel'),'String','$Cart~Position~$(meters)','interpreter','latex')
% set(get(AX(2),'Ylabel'),'String','$Pendulum~Angle~$(radians)','interpreter','latex')
% xlabel('$Time~$(seconds)','interpreter','latex')
% 
% title('Step Response with Observer-Based State-Feedback Control')

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