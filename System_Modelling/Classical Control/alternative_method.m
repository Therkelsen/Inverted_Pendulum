clc, clear, close all

M = 0.5; % mass of pend [kg]
m = 0.084;  % mass of cart [kg]
b = 5; % viscous dampening of cart [N/(m/s)]
g = 9.82; % acceleration of gravity [m/s^2]
l = 0.35; % length of pend [m]
I = (1.0/3.0)*m*l^2; % Moment of inertia for rod

q = (M+m)*(I+m*l^2)-(m*l)^2;
s = tf('s');

% Transfer functions for the system are modelled as per
% https://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling

tf_cart = (((I+m*l^2)/q)*s^2 - (m*g*l/q))/(s^4 + (b*(I + m*l^2))*s^3/q - ((M + m)*m*g*l)*s^2/q - b*m*g*l*s/q)
tf_pend = (m*l*s/q)/(s^3 + (b*(I + m*l^2))*s^2/q - ((M + m)*m*g*l)*s/q - b*m*g*l/q)
y = m*l/q
s = tf('s');

tf_sys = [tf_cart ; tf_pend];

inputs = {'u'};
outputs = {'x'; 'phi'};

set(tf_sys,'InputName',inputs)
set(tf_sys,'OutputName',outputs)

tf_sys

% System is shown on state space form
ss_sys = ss(tf_sys)

% The system is analysed as per
% https://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemAnalysis

t=0:0.01:1;
figure
impulse(tf_sys,t);
title('Open-Loop Impulse Response')

figure
pzmap(tf_cart,tf_pend)
grid on
[zeros_c, poles_p] = zpkdata(tf_cart,'v')
[zeros_p, poles_p] = zpkdata(tf_pend,'v')

t = 0:0.05:10;
u = ones(size(t));
[y,t] = lsim(tf_sys,u,t);
figure
plot(t,y)
title('Open-Loop Step Response')
axis([0 3 0 50])
legend('x','phi')

step_info = lsiminfo(y,t);
cart_info = step_info(1)
pend_info = step_info(2)

%% PID Tuning
close all
[C_pi,info] = pidtune(tf_pend,'PIDF')

tf_pend_pi = feedback(C_pi*tf_pend, 1)

stepinfo(tf_pend_pi)

figure()
step(tf_pend_pi)

figure()
impulse(tf_pend_pi)

figure()
margin(C_pi*tf_pend)

%%
% % Now, root locus is used to find the proper controller as per
% % https://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=ControlRootLocus
% 
% figure
% rlocus(tf_pend)
% title('Root Locus of Plant (under Proportional Control)')
% 
% C = 1/s;
% figure
% rlocus(C*tf_pend)
% title('Root Locus with Integral Control')
% 
% zeros = zero(C*tf_pend)
% poles = pole(C*tf_pend)
% 
% % PID Controller is made
% z = [-3 -4];
% p = 0;
% k = 1;
% C = zpk(z,p,k);
% 
% rlocus(C*tf_pend)
% title('Root Locus with PID Controller')
% [k,poles] = rlocfind(C*tf_pend)
% 
% K = k;
% T = feedback(tf_pend,K*C);
% figure
% impulse(T)
% title('Impulse Disturbance Response of Pendulum Angle under PID Control');
% 
% T2 = feedback(1,tf_pend*C)*tf_cart;
% t = 0:0.01:8.5;
% figure
% impulse(T2, t);
% title('Impulse Disturbance Response of Cart Position under PID Control');
% 
% step(T,t)