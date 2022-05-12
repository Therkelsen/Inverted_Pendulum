clc, clear, close all
syms mc mp l theta dtheta ddtheta x dx ddx F bc bp Jp g

g_ = 9.82    %Earth gravitational force

mc_ = 0.5    %Mass of cart
bc_ = 5      %Friction coeficient of cart

mp_ = 0.084  %Mass of pendulum
bp_ = 0.0012 %Friction coeficcient of pendulum
l_ = 0.35    %length of pendulum
Jp_ = (1/3)*mp_*l_^2

%Ours (probably wrong)
%cart = (mc+mp)*ddx + mp*l*cos(theta)*ddtheta - mp*l*sin(theta)*dtheta^2 == F - bc*dx
%pend = (mp*l + Jp)*ddtheta + mp*l*ddx*cos(theta)+mp*g*l*sin(theta) == 0

cart = (mc+mp)*ddx-mp*l*cos(theta)*ddtheta+mp*l*sin(theta)*dtheta^2 == F - bc*dx
pend = (mp*l + Jp)*ddtheta-mp*l*ddx*cos(theta)-mp*g*l*sin(theta) == -bp*dtheta



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

Al = subs(An, {theta dtheta x dx F},{0 0 0 0 0})
Bl = subs(Bn, {theta dtheta x dx F},{0 0 0 0 0})

%Insert Values
A = subs(Al, {Jp mc mp l bc bp g}, {Jp_ mc_ mp_ l_ bc_ bp_ g_});
B = subs(Bl, {Jp mc mp l bc bp g}, {Jp_ mc_ mp_ l_ bc_ bp_ g_});

vpa(A)
vpa(B)

%Values for matrix [cart pos, cart velocity, pend pos, pend velocity]
%velocity
C = [1 0 0 0; 0 0 1 0] %Cart speed
D = 0;


% Statespace object
sys_ss = ss(double(A),double(B),double(C),D);

% Transferfunction from statespace
sys_tf = tf(sys_ss)

% Individual transfer functions
[num, denum] = tfdata(sys_tf);

sys_tf_cart = tf(num(1),denum(1));
sys_tf_pend = tf(num(2),denum(2));



