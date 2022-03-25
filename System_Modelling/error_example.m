clc; clear; close all;

% mPend = 0.084; % mass of pend [kg]
% M = 0.5; % mass of cart [kg]
% mRod = 0.082;  % mass of rod [kg]
% m = mRod + mPend;
% L = 0.35; % length of pend [m]
% g = -9.82; % acceleration of gravity [m/s^2]
% bCart = 5; % viscous dampening of cart [N/(m/s)]
% bPend = 0.0012; % viscious dampening coeff for pend[N/(m/s)]
% I = (1.0/3.0)*m*L^2; % Moment of inertia for rod

syms x(t) theta(t) M m b t l

% Sum of horizontal forces on pendulum
F_pend_hor = m*diff(diff(x,t),t) + m*l*diff(diff(theta,t),t) + cos(theta) - m*l*diff(theta,t)^2*sin(theta)

% Sum of forces of the cart in the horizontal direction, which needs to
% have the motion of the pendulum added to the end
F_cart_hor = M*diff(diff(x,t),t)

F_sys = F_pend_hor + F_cart_hor